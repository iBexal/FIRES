######################## SETUP ########################
import os

from fires_base import *
from order_tracing import *

current_folder = os.path.dirname(os.path.abspath(__file__))

#######################################################

# EXAMPLE USAGE OF FIRES
data_folder = os.path.join(current_folder, 'test_data', 'correct') # Give path of files to check. check_files() will look for all .fit files in this folder.
check_files(data_folder, log_obs_type=True, log_new_coords=True) #This will check all files in the folder, save new .fit files in a subfolder of data_folder, and save a log file in the same subfolder.
# END EXAMPLE