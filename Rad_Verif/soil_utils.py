#!/bin/python

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import netCDF4
from cdo import * #This is bad style, unfortunately pyCDO was implemented in a
# way that would make anything else too cumbersome.

class Verify:
    '''This class contains tools for data I/O, remapping and score calcualtions.
    initialize with model and obs paths'''
    __author__ = "Daniel Regenass, daniel.regenass@env.ethz.ch"


    def __init__(self, model_path, obs_path, **kwargs):

        self.model_path = model_path
        self.obs_path = obs_path

        #initialize CDO operators
        self.cdo = Cdo()

        #choose variables if specified
        if "model_names" in kwargs and "obs_names" in kwargs:
            self.model_names = kwargs.get("model_names")
            self.obs_names = kwargs.get("obs_names")
            print("Choosing variables %s" %self.model_names)
            self.cdo.selname(self.model_names, input=self.model_path, output="tmp1.nc")
            self.cdo.selname(self.obs_names, input=self.obs_path, output="tmp2.nc")
            #change file path for next step
            self.model_path = "tmp1.nc"
            self.obs_path = "tmp2.nc"
        else:
            print("Reading in full files (all variables) as dataset")

        #remap with cdo if field to be remapped is specified
        if "remap_type" in kwargs:
            self.remap_type = kwargs.get("remap_type")
            if "grid" in kwargs:
                self.grid = kwargs.get("grid")
            else:
                self.grid = None
            print("Remapping %s"  %self.remap_type)
            self.__remap(grid=self.grid)

        #read data
        self.model_data = netCDF4.Dataset(model_path)
        self.obs_data = netCDF4.Dataset(obs_path)

        #set reference data field for attributes,etc.
        if self.remap_type == "mod":
            self.ref_data = self.model_data
            self.ref_names = self.model_names
        else:
            self.ref_data = self.obs_data
            self.ref_names = self.obs_names

        # automatic extraction of fields, if only one variable is given
        nondim_vars_in_input = self.model_names.split(",")
        if len(nondim_vars_in_input) == 1:
            self.extract_fields(self.model_names, self.obs_names)
        else:
            raise ValueError("More then one input field is currently not\
supported! Aborting now ..")




    def __remap(self, **kwargs):

        if self.remap_type == "obs":
            field_to_remap = self.obs_path
            reference_field = self.model_path
            remapped_field = "regridded_" + self.obs_path
        elif self.remap_type == "mod":
            field_to_remap = self.model_path
            reference_field = self.obs_path
            remapped_field = "regridded_" + self.model_path

        # check if grid is provided, else extract it from reference field.
        if "grid" in kwargs:
            grid = kwargs.get("grid")
        else:
            grid = "extracted_grid.txt"
            print(grid)
            # extract output (stdout) from griddes
            with open(grid, "w") as sys.stdout:
                self.cdo.griddes(input=reference_field)

        # bilinear remapping
        self.cdo.remapbil(grid, input=field_to_remap, output=remapped_field)




    def extract_fields(self, fieldname_mod, fieldname_obs):
        '''Extract fields (variables) by name, one for model data and one for
        obs data. Creates the tuple self.fields'''

        print(fieldname_mod)

        mod_field = self.model_data.variables[fieldname_mod][:]
        obs_field = self.obs_data.variables[fieldname_obs][:]

        self.out_dims = self.ref_data.variables[self.ref_names].dimensions
        self.out_datatype = self.ref_data.variables[self.ref_names].datatype
        self.out_attr = \
        {k: self.ref_data.variables[self.ref_names].getncattr(k) \
          for k in self.ref_data.variables[self.ref_names].ncattrs()}

        self.__fields = (mod_field, obs_field)


    def restrict_area(self, )



    def calculate_bias(self):
        '''Calculates bias (mean error)'''
        self.bias = self.__fields[0] - self.__fields[1]



    def calculate_rmse(self):
        '''Calcualtes RMSE'''
        self.rmse = np.sqrt((self.__fields[0] - self.__fields[1])**2)



    def calculate_sd(self):
        '''Calculate Error SD (prefferably from RMSE and ME)'''
        self.sd = self.rmse - self.bias


    def calculate_scores(self):
        '''ME, RMSE and SD'''
        self.calculate_bias()
        self.calculate_rmse()
        self.calculate_sd()


    def scores_to_file(self, scorenames, filename):
        '''Write chosen scores to netCDF4 file'''

        print("Creating output file.")
        scorefile = netCDF4.Dataset(filename, "w")

        #Copy dimensions
        print("Copying dimensions from original file.")
        for dim_name, the_dim in self.ref_data.dimensions.iteritems():
            print dim_name, len(the_dim)
            scorefile.createDimension(dim_name, len(the_dim) \
              if not the_dim.isunlimited() else None)
            # also as variables if they exist
            if dim_name in self.ref_data.variables:
                dim_var = scorefile.createVariable(dim_name, \
                  self.ref_data.variables[dim_name].datatype, \
                  self.ref_data.variables[dim_name].dimensions)
                dim_attr =  \
                  {k: self.ref_data.variables[dim_name].getncattr(k) \
                  for k in self.ref_data.variables[dim_name].ncattrs()}
                dim_var.setncatts(dim_attr)
                dim_var[:] = self.ref_data.variables[dim_name][:]


        print("Create variables with attributes")
        # Copy variable attributes

        score_dict = {"ME":self.bias, "RMSE":self.rmse,"SD":self.sd}

        for scorename in scorenames:

            print(scorename)

            if not scorename in score_dict:
                break

            out_var = scorefile.createVariable(scorename, self.out_datatype,\
              self.out_dims)
            out_var.setncatts(self.out_attr)#\
            #      for k in self.mod_field.ncattrs()})
            out_var[:] = score_dict[scorename]

        scorefile.close()
        print("Done writing data.")
