#Assigns Base Url

base_url <- "http://navigator.oceansdata.ca/api/v1.0/subset/"


#Assignes Query Variables
dataset_name= '"riops_monthly"'
max_range= '"64.11624382643063,-38.932616867971234"'
min_range= '"49.937662657046786,-61.015823044729046"'
output_format= '"NETCDF4"'
should_zip= 0
time= '"13,25"'
user_grid= 0
variables= '"vozocrtx,vomecrty"'

#Assembles Query Using the Assigned Values
query = sprintf('{"dataset_name": %s, "max_range": %s, "min_range": %s, "output_format": %s, "should_zip": %s, "time": %s, "user_grid": %s, "variables": %s}', dataset_name, max_range, min_range, output_format, should_zip, time, user_grid, variables)

#Request and Save Image
full_url <- paste0(base_url, "?query=", URLencode(query, reserved=TRUE))
#Format time to be used in file name
time_ = gsub(':.*','', time)
time_ = gsub('\"',"",time_)

filename = paste0("script_output_", time_, ".nc")
download.file(full_url, filename, extra = "curl", mode = "wb")
