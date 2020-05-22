# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 22:34:36 2018

@author: kkc29
"""

from __future__ import print_function

from PYME.recipes.base import ModuleBase, register_module, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrList, File

import numpy as np
import pandas as pd
from PYME.IO import tabular, MetaDataHandler
from PYME.IO.image import ImageStack
import os
from collections import OrderedDict
from matplotlib.pyplot import *
from scipy import optimize, interpolate
from functools import partial

class CalculateFRCBase(ModuleBase):
    """
    Base class. Refer to derived classes for docstrings.
    """
    pre_filter = Enum(['Tukey_1/8', None])
    frc_smoothing_func = Enum(['Cubic Spline', 'Sigmoid', None])
    multiprocessing =  Bool(True)
    plot_graphs = Bool()
    
#    output_fft_image_a = Output('FRC_fft_image_a')
#    output_fft_image_b = Output('FRC_fft_image_b')
    output_fft_images_cc = Output('FRC_fft_images_cc')
    output_frc_dict = Output('FRC_dict')
    output_frc_plot = Output('FRC_plot')
    
    def execute(self):
        raise Exception("Base class not fully implemented")
    
    def preprocess_images(self, image_pair):
        # pad images to square shape
#        image_pair = self.pad_images_to_equal_dims(image_pair)
        
        dims_length = np.stack([im.shape for im in image_pair], 0)
        assert np.all([np.all(dims_length[:, i] == dims_length[0, i])for i in xrange(dims_length.shape[1])]), "Images not the same dimension."
        
        # apply filtering to zero near edge of images
        if self.pre_filter == 'Tukey_1/8':
            image_pair = self.filter_images_tukey(image_pair, 1./8)
        elif self.pre_filter == None:
            pass
        else:
            raise Exception()
            
        return image_pair
        
#    def pad_images_to_equal_dims(self, images):
#        dims_length = np.stack([im.shape for im in images], 0)
#        assert np.all([np.all(dims_length[:, i] == dims_length[0, i])for i in xrange(dims_length.shape[1])]), "Images not the same dimension."
#        
#        return images
#        
#        dims_length = dims_length[0, :]        
#        max_dim = dims_length.max()
#        
#        padding = np.empty((dims_length.shape[0], 2), dtype=np.int)
#        for dim in xrange(dims_length.shape[0]):
#            total_padding = max_dim - dims_length[dim]
#            padding[dim] = [total_padding // 2, total_padding - total_padding //2]
#            
#        results = list()
#        for im in images:
#            results.append(np.pad(im, padding, mode='constant', constant_values=0))
#        
#        return results

    def filter_images_tukey(self, images, alpha):
        from scipy.signal import tukey
        
#        window = tukey(images[0].shape[0], alpha=alpha)
#        window_nd = np.prod(np.stack(np.meshgrid(*(window,)*images[0].ndim)), axis=0)
        windows = [tukey(images[0].shape[i], alpha=alpha) for i in xrange(images[0].ndim)]
        window_nd = np.prod(np.stack(np.meshgrid(*windows, indexing='ij')), axis=0)
        
#        for im in images:
#            im *= window_nd
        
        return [images[0]*window_nd, images[1]*window_nd]

    def calculate_FRC_from_images(self, image_pair, mdh):
        ft_images = list()
        if self.multiprocessing:
            results = list()
            for im in image_pair:
                results.append(self._pool.apply_async(np.fft.fftn, (im,)))            
            for res in results:
                ft_images.append(res.get())
            del results
        else:
            for im in image_pair:
                ft_images.append(np.fft.fftn(im))
        
#        im_fft_freq = np.fft.fftfreq(image_pair[0].shape[0], self._pixel_size_in_nm)
#        im_R = np.sqrt(im_fft_freq[:, None]**2 + im_fft_freq[None, :]**2)
        im_fft_freqs = [np.fft.fftfreq(image_pair[0].shape[i], self._pixel_size_in_nm[i]) for i in xrange(image_pair[0].ndim)]
        im_R = np.linalg.norm(np.stack(np.meshgrid(*im_fft_freqs, indexing='ij')), axis=0)

        im1_fft_power = np.multiply(ft_images[0], np.conj(ft_images[0]))
        im2_fft_power = np.multiply(ft_images[1], np.conj(ft_images[1]))        
        im12_fft_power = np.multiply(ft_images[0], np.conj(ft_images[1]))                
        
##        fft_ims = ImageStack(data=np.stack([np.fft.fftshift(im1_fft_power),
##                                            np.fft.fftshift(im2_fft_power),
##                                            np.fft.fftshift(im12_fft_power)], axis=-1), mdh=mdh)
##        self._namespace[self.output_fft_images] = fft_ims
#        self._namespace[self.output_fft_image_a] = ImageStack(data=np.fft.fftshift(im1_fft_power), titleStub="ImageA_FFT")
#        self._namespace[self.output_fft_image_b] = ImageStack(data=np.fft.fftshift(im2_fft_power), titleStub="ImageB_FFT")
#        self._namespace[self.output_fft_images_cc] = ImageStack(data=np.fft.fftshift(im12_fft_power), titleStub="ImageA_Image_B_FFT_CC")
        
        self._namespace[self.output_fft_images_cc] = ImageStack(data=np.stack([np.atleast_3d(np.fft.fftshift(im1_fft_power)),
               np.atleast_3d(np.fft.fftshift(im2_fft_power)),
               np.atleast_3d(np.fft.fftshift(im12_fft_power))], 3), titleStub="ImageA_Image_FFT_CC")
            
        if self.plot_graphs:
            from PYME.DSView.dsviewer import ViewIm3D, View3D
#            ViewIm3D(self._namespace[self.output_fft_image_a])
#            ViewIm3D(self._namespace[self.output_fft_image_b])
            ViewIm3D(self._namespace[self.output_fft_images_cc])
#            View3D(np.fft.fftshift(im_R))
            
        
        im1_fft_flat_res = CalculateFRCBase.BinData(im_R.flatten(), im1_fft_power.flatten(), statistic='mean', bins=201)
        im2_fft_flat_res = CalculateFRCBase.BinData(im_R.flatten(), im2_fft_power.flatten(), statistic='mean', bins=201)
        im12_fft_flat_res = CalculateFRCBase.BinData(im_R.flatten(), im12_fft_power.flatten(), statistic='mean', bins=201)
        
        corr = np.real(im12_fft_flat_res.statistic) / np.sqrt(np.abs(im1_fft_flat_res.statistic*im2_fft_flat_res.statistic))
        
        smoothed_frc = self.smooth_frc(im12_fft_flat_res.bin_edges[:-1], corr)
        
        res = self.calculate_threshold(im12_fft_flat_res.bin_edges[:-1], corr, smoothed_frc, im12_fft_flat_res.counts)
        
        return res
    
    def smooth_frc(self, freq, corr):
        if self.frc_smoothing_func is None:
            interp_frc = interpolate.interp1d(freq, corr, kind='next', )
            return interp_frc
        
        elif self.frc_smoothing_func == "Sigmoid":
            func = CalculateFRCBase.Sigmoid
            fit_res = optimize.minimize(lambda a, x: np.sum(np.square(func(a[0], a[1], x)-corr)), [1, freq[len(freq)/2]], args=(freq), method='Nelder-Mead')

            return partial(func, n=fit_res.x[0], c=fit_res.x[1])
        elif self.frc_smoothing_func == "Cubic Spline":
            # smoothed so that average deviation loss is less than 0.2% of original. Somewhat arbitrary but probably not totally unreasonable since FRC is bounded 0 to 1.
            interp_frc = interpolate.UnivariateSpline(freq, corr, k=3, s=len(freq)*(0.002*np.var(corr)))

            return interp_frc            
    
    def calculate_threshold(self, freq, corr, corr_func, counts):
        res = dict()
        
        fsc_0143 = optimize.minimize(lambda x: np.square(corr_func(x=x)-0.143), freq[np.argmax(corr_func(x=freq)-0.143 < 0)], method='Nelder-Mead')
        res['frc 1/7'] = 1./fsc_0143.x[0]
        
        sigma = 1.0 / np.sqrt(counts*0.5)
        sigma_spl = interpolate.UnivariateSpline(freq, sigma, k=3, s=0)
        fsc_3sigma = optimize.minimize(lambda x: np.square(corr_func(x=x)-3.*sigma_spl(x)), freq[np.argmax(corr_func(x=freq)-3.*sigma_spl(freq) < 0)], method='Nelder-Mead')
        res['frc 3 sigma'] = 1./fsc_3sigma.x[0]        
        
        # van Heel and Schatz, 2005, Fourier shell correlation threshold criteria    
        half_bit = (0.2071 + 1.9102 / np.sqrt(counts)) / (1.2071 + 0.9102 / np.sqrt(counts))
        half_bit_spl = interpolate.UnivariateSpline(freq, half_bit, k=3, s=0)
        fsc_half_bit = optimize.minimize(lambda x: np.square(corr_func(x=x)-half_bit_spl(x)), freq[np.argmax(corr_func(x=freq)-half_bit_spl(freq) < 0)], method='Nelder-Mead')
        res['frc half bit'] = 1./fsc_half_bit.x[0]
        
#        fsc_max = np.max([fsc_0143.x[0], fsc_2sigma.x[0], fsc_3sigma.x[0], fsc_5sigma.x[0], fsc_half_bit.x[0]])
#        axes[1].set_xlim(0, np.min([2*fsc_max, im12_fft_flat_res.bin_edges[-1]]))
        
        if not self.plot_graphs:
            ioff()
            
        frc_text = ""
        fig, axes = subplots(1,2,figsize=(10,4))
        axes[0].plot(freq, corr)
        axes[0].plot(freq, corr_func(x=freq))
        
        axes[0].axhline(0.143, ls='--', color='red')
        axes[0].axvline(fsc_0143.x[0], ls='--', color='red', label='1/7')
        frc_text += "\nFRC 1/7:  {:.2f} nm".format(1./fsc_0143.x[0])
 
        axes[0].plot(freq, 3*sigma_spl(freq), ls='--', color='pink')
        axes[0].axvline(fsc_3sigma.x[0], ls='--', color='pink', label='3 sigma')
        frc_text += "\nFRC 3 sigma:  {:.2f} nm".format(1./fsc_3sigma.x[0])
 
        axes[0].plot(freq, half_bit_spl(freq), ls='--', color='purple')
        axes[0].axvline(fsc_half_bit.x[0], ls='--', color='purple', label='1/2 bit')

        frc_text += "\n1/2 bit:  {:.2f} nm".format(1./fsc_half_bit.x[0])
        
        axes[0].legend()            
#        axes[0].set_ylim(None, 1.1)
        
        x_ticklocs = axes[0].get_xticks()
        axes[0].set_xticklabels(["{:.1f}".format(1./i) for i in x_ticklocs])
        axes[0].set_ylabel("FRC")
        axes[0].set_xlabel("Resol (nm)")
        
        axes[1].text(0.5, 0.5, frc_text, horizontalalignment='center', verticalalignment='center', transform=axes[1].transAxes)
        axes[1].set_axis_off()
        
        if self.plot_graphs:
            fig.show()
        else:
            ion()
        
        self._namespace[self.output_frc_plot] = fig
        
        return res
    
    @staticmethod
    def BinData(indexes, data, statistic='mean', bins=10):
        # Calculates binned statistics. Supports complex number.
        if statistic == 'mean':
            func = np.mean
        elif statistic == 'sum':
            func = np.sum
        
        class Result(object):
            statistic = None
            bin_edges = None
            counts = None
        
        bins = np.linspace(indexes.min(), indexes.max(), bins)
        binned = np.zeros(len(bins)-1, dtype=data.dtype)
        counts = np.zeros(len(bins)-1, dtype=np.int)
        
        indexes_sort_arg = np.argsort(indexes.flatten())
        indexes_sorted = indexes.flatten()[indexes_sort_arg]
        data_sorted = data.flatten()[indexes_sort_arg]
        edges_indexes = np.searchsorted(indexes_sorted, bins)
        
        for i in xrange(bins.shape[0]-1):
            values = data_sorted[edges_indexes[i]:edges_indexes[i+1]]
            binned[i] = func(values)
            counts[i] = len(values)
            
        res = Result()
        res.statistic = binned
        res.bin_edges = bins
        res.counts = counts
        return res
    
    @staticmethod
    def Sigmoid(n, c, x):
        res = 1 - 1 / (1 + np.exp(n*(-x+c)))
        return res
    
@register_module('FSCFromImages')
class CalculateFRCFromImages(CalculateFRCBase):
    """
    Take a pair of images and calculates the fourier shell/ring correlation (FSC / FRC).
    
    Inputs
    ------
    input_image_a : ImageStack
        First of two images.
        
    Outputs
    -------
    output_fft_image_a : ImageStack
        Fast Fourier transform of the first image.
    output_fft_image_b : ImageStack
        Fast Fourier transform of the second image.
    output_fft_images_cc : ImageStack
        Fast Fourier transform cross-correlation.
    output_frc_dict : dict
        FSC/FRC results.
    output_frc_plot : Plot
        Output plot of the FSC / FRC curve.
    
    Parameters
    ----------
    image_b_path : File
        File path of the second of the two images.
    c_channel : int
        Color channel of the images to use.
    image_a_z : int
        Ignored unless flatten_z is True. In which case either select the z plane to use (>=0) or performs a maximum project (<0) for the first image.
    image_b_z : int
        Ignored unless flatten_z is True. In which case either select the z plane to use (>=0) or performs a maximum project (<0) for the second image.
    flatten_z : Bool
        If enabled ignores z information and only performs a FRC.
    pre_filter : string
        Methods to filter the images prior to Fourier transform.
    frc_smoothing_func : string
        Methods to smooth the FSC / FRC curve.
    multiprocessing : Bool
        Enables multiprocessing.
    plot_graphs : Bool
        Show graphs.
    
    """
    
    input_image_a = Input('input')
#    image_a_dim = Int(2)
#    image_a_index = Int(0)
#    image_b_dim = Int(2)
#    image_b_index = Int(1)
    image_b_path = File()
    c_channel = Int(0)
    flatten_z = Bool(True)
    image_a_z = Int(-1)
    image_b_z = Int(-1)
    
    def execute(self, namespace):
        self._namespace = namespace
        import multiprocessing
        from PYME.util import mProfile
        
        mProfile.profileOn(["frc.py"])
        
        if self.multiprocessing:
            proccess_count = np.clip(2, 1, multiprocessing.cpu_count()-1)
            self._pool = multiprocessing.Pool(processes=proccess_count)        
       
#        image_pair = self.generate_image_pair(mapped_pipeline)
#        ims = namespace[self.input_images]
        image_a = namespace[self.input_image_a]
        image_b = ImageStack(filename=self.image_b_path)
        
        self._pixel_size_in_nm = np.zeros(3, dtype=np.float)
        self._pixel_size_in_nm[0] = image_a.mdh.voxelsize.x
        self._pixel_size_in_nm[1] = image_a.mdh.voxelsize.y
        try:
            self._pixel_size_in_nm[2] = image_a.mdh.voxelsize.z
        except:
            pass
        if image_a.mdh.voxelsize.units == 'um':
            self._pixel_size_in_nm *= 1.E3
            print(self._pixel_size_in_nm)
        
#        image_indices = [[self.image_a_dim, self.image_a_index], [self.image_b_dim, self.image_b_index]]
#        image_slices = list()
#        for i in xrange(2):
#            slices = [slice(None, None), slice(None, None)]
#            for j in xrange(2, image_indices[i][0]+1):
#                if j == image_indices[i][0]:
#                    slices.append(slice(image_indices[i][1], image_indices[i][1]+1))
#                else:
#                    slices.append(slice(None, None))
#            image_slices.append(slices)
#        
#        image_pair = [ims.data[image_slices[0]].squeeze(), ims.data[image_slices[1]].squeeze()]
        image_a_data = image_a.data[:,:,:,self.c_channel].squeeze()
        image_b_data = image_b.data[:,:,:,self.c_channel].squeeze()
        if self.flatten_z:
            print("2D mode. Slice if z index >= 0 otherwise max projection")
            if self.image_a_z >= 0:
                image_a_data = image_a_data[:,:,self.image_a_z]
            else:
                image_a_data = image_a_data.max(2)
            if self.image_b_z >= 0:
                image_b_data = image_b_data[:,:,self.image_b_z]
            else:
                image_b_data = image_b_data.max(2)
#            print(np.allclose(image_a_data, image_b_data))
        image_pair = [image_a_data, image_b_data]
#        print(image_pair[0].shape)
        
        image_pair = self.preprocess_images(image_pair)            
       
        frc_res = self.calculate_FRC_from_images(image_pair, None)
        
        namespace[self.output_frc_dict] = frc_res
        
        if self.multiprocessing:
            self._pool.close()
            self._pool.join()
        
        mProfile.profileOff()
        mProfile.report()
    

@register_module('FSCFromLocs')
class CalculateFRCFromLocs(CalculateFRCBase):
    """
    Generates a pair of images from localization data and calculates the fourier shell/ring correlation (FSC / FRC).
    
    Inputs
    ------
    inputName : TabularBase
        Localization data.
        
    Outputs
    -------
    outputName : TabularBase
        Localization data labeled with how it was divided (FRC_group).
    output_images : ImageStack
        Pair of 2D or 3D histogram rendered images.
    output_fft_image_a : ImageStack
        Fast Fourier transform of the first image.
    output_fft_image_b : ImageStack
        Fast Fourier transform of the second image.
    output_fft_images_cc : ImageStack
        Fast Fourier transform cross-correlation.
    output_frc_dict : dict
        FSC/FRC results.
    output_frc_plot : Plot
        Output plot of the FSC / FRC curve.
    
    Parameters
    ----------
    split_method : string
        Different methods of dividing data into halves.
    pixel_size_in_nm : int
        Pixel size used for rendering the images.
    flatten_z : Bool
        If enabled ignores z information and only performs a FRC.
    pre_filter : string
        Methods to filter the images prior to Fourier transform.
    frc_smoothing_func : string
        Methods to smooth the FSC / FRC curve.
    multiprocessing : Bool
        Enables multiprocessing.
    plot_graphs : Bool
        Show graphs.
    
    """
    inputName = Input('Localizations')
    split_method = Enum(['halves_random', 'halves_time', 'halves_100_time_chunk', 'halves_10_time_chunk', 'fixed_time', 'fixed_10_time_chunk'])    
    pixel_size_in_nm = Int(5)
#    pre_filter = Enum(['Tukey_1/8', None])
#    frc_smoothing_func = Enum(['Cubic Spline', 'Sigmoid', None])
#    plot_graphs = Bool()
#    multiprocessing =  Bool(True)
    flatten_z = Bool(True)
    
    outputName = Output('FRC_ds')
    output_images = Output('FRC_images')
#    output_fft_images = Output('FRC_fft_images')
#    output_frc_dict = Output('FRC_dict')
#    output_frc_plot = Output('FRC_plot')
    
    def execute(self, namespace):
        self._namespace = namespace
        import multiprocessing
        from PYME.util import mProfile
        
        mProfile.profileOn(["frc.py"])
        
        if self.multiprocessing:
            proccess_count = np.clip(2, 1, multiprocessing.cpu_count()-1)
            self._pool = multiprocessing.Pool(processes=proccess_count)
        
        pipeline = namespace[self.inputName]
        mapped_pipeline = tabular.mappingFilter(pipeline)
        self._pixel_size_in_nm = self.pixel_size_in_nm * np.ones(3, dtype=np.float)
        
        image_pair = self.generate_image_pair(mapped_pipeline)
        
        image_pair = self.preprocess_images(image_pair)
            
        # Should use DensityMapping recipe eventually when it is ready.
        mdh = MetaDataHandler.NestedClassMDHandler()
        mdh['Rendering.Method'] = "np.histogramdd"
        if 'imageID' in pipeline.mdh.getEntryNames():
            mdh['Rendering.SourceImageID'] = pipeline.mdh['imageID']
        try:
            mdh['Rendering.SourceFilename'] = pipeline.resultsSource.h5f.filename
        except:
            pass        
        mdh.Source = MetaDataHandler.NestedClassMDHandler(pipeline.mdh)        
        mdh['Rendering.NEventsRendered'] = [image_pair[0].sum(), image_pair[1].sum()]
        mdh['voxelsize.units'] = 'um'
        mdh['voxelsize.x'] = self.pixel_size_in_nm * 1E-3
        mdh['voxelsize.y'] = self.pixel_size_in_nm * 1E-3
        
        ims = ImageStack(data=np.stack(image_pair, axis=-1), mdh=mdh)
        namespace[self.output_images] = ims
        
        if self.plot_graphs:
            from PYME.DSView.dsviewer import ViewIm3D
            ViewIm3D(ims)
        
        frc_res = self.calculate_FRC_from_images(image_pair, pipeline.mdh)
        
#        smoothed_frc = self.SmoothFRC(frc_freq, frc_corr)
#        
#        self.CalculateThreshold(frc_freq, frc_corr, smoothed_frc)
        
        namespace[self.output_frc_dict] = frc_res
        
        if self.multiprocessing:
            self._pool.close()
            self._pool.join()
        
        mProfile.profileOff()
        mProfile.report()
        
    def generate_image_pair(self, mapped_pipeline):
        # Split localizations into 2 sets
        mask = np.zeros(mapped_pipeline['t'].shape, dtype=np.bool)        
        if self.split_method == 'halves_time':
            sort_arg = np.argsort(mapped_pipeline['t'])
            mask[sort_arg[:len(sort_arg)/2]] = 1
        elif self.split_method == 'halves_random':
            mask[:len(mask)/2] = 1
            np.random.shuffle(mask)            
        elif self.split_method == 'halves_100_time_chunk':
            sort_arg = np.argsort(mapped_pipeline['t'])
            chunksize = mask.shape[0] / 100.0
            for i in range(50):
                mask[sort_arg[int(np.round(i*2*chunksize)):int(np.round((i*2+1)*chunksize))]] = 1
        elif self.split_method == 'halves_10_time_chunk':
            sort_arg = np.argsort(mapped_pipeline['t'])
            chunksize = mask.shape[0] * 0.1
            for i in range(5):
                mask[sort_arg[int(np.round(i*2*chunksize)):int(np.round((i*2+1)*chunksize))]] = 1
        elif self.split_method == 'fixed_time':
            time_cutoff = (mapped_pipeline['t'].ptp() + 1) // 2 + mapped_pipeline['t'].min()
            mask[mapped_pipeline['t'] < time_cutoff] = 1
        elif self.split_method == 'fixed_10_time_chunk':
            time_cutoffs = np.linspace(mapped_pipeline['t'].min(), mapped_pipeline['t'].max(), 11, dtype=np.float)
            for i in xrange((time_cutoffs.shape[0]-1)//2):
                mask[(mapped_pipeline['t'] > time_cutoffs[i*2]) & (mapped_pipeline['t'] < time_cutoffs[i*2+1])] = 1
        
        if self.split_method.startswith('halves_'):
            assert np.abs(mask.sum() - (~mask).sum()) <= 1, "datasets uneven, {} vs {}".format(mask.sum(), (~mask).sum())
        else:
            print("variable counts between images, {} vs {}".format(mask.sum(), (~mask).sum()))
        
        mapped_pipeline.addColumn("FRC_group", mask)
        self._namespace[self.outputName] = mapped_pipeline        
        
        dims = ['x', 'y', 'z']
        if self.flatten_z:
            dims.remove('z')
            
        # Simple hist2d binning
        bins = list()
        data = list()
        for i, dim in enumerate(dims):
            bins.append(np.arange(np.floor(mapped_pipeline[dim].min()).astype(int), np.ceil(mapped_pipeline[dim].max()).astype(int) + self.pixel_size_in_nm, self.pixel_size_in_nm))
            
            data.append(mapped_pipeline[dim])
        
        data = np.stack(data, axis=1)
        print(bins)
        
        if self.multiprocessing:
            results = list()
            results.append(self._pool.apply_async(np.histogramdd, (data[mask==0],), kwds= {"bins":bins}))
            results.append(self._pool.apply_async(np.histogramdd, (data[mask==1],), kwds= {"bins":bins}))
            
            image_a = results[0].get()[0]
            image_b = results[1].get()[0]                
        else:
            image_a = np.histogramdd(data[mask==0], bins=bins)[0]
            image_b = np.histogramdd(data[mask==1], bins=bins)[0]
            
        return image_a, image_b