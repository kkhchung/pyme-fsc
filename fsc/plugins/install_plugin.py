from PYME import config
import os
import sys
import pkgutil
import recipes, dsviewer, visgui # the 3 types of plugins

def main():
#    this_dir = os.path.dirname(__file__)

    if len(sys.argv) > 1 and sys.argv[1] == 'dist':
        target_dir = config.dist_config_directory
    else:  # no argument provided or is not 'dist', default to user config directory
        target_dir = config.user_config_dir    

    target_filename =  recipes.__name__.split(".")[0] + ".txt"
    for plugin in [recipes, dsviewer, visgui]:
        target_sub_dir = plugin.__name__.split(".")[-1]
        target_path = os.path.join(target_dir, "plugins", target_sub_dir, target_filename)
        modules_list = create_module_list(plugin)
        if len(modules_list) > 0:
            print("writing addons to file: {}".format(target_path))
            for val in modules_list:
                print val
            with open(target_path, 'w') as f:
                f.writelines(modules_list)
        
def create_module_list(plugin):
    modules = list()    
    for _, name, _ in pkgutil.iter_modules([os.path.dirname(plugin.__file__)]):
        modules.append(".".join([plugin.__name__, name])+'\n')
    return modules

if __name__ == '__main__':
    main()
#    create_module_list()