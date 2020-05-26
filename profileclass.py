from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
import numpy as np


class Profile:

    def __init__(self, profile_mat):
        self.original_positions = profile_mat[:,0]
        self.original_values = profile_mat[:,1]
        self.position = profile_mat[:,0]
        self.profilevalues = profile_mat[:,1]
        self.xscale = 1.0
        self.yscale = 1.0
        self.xoffset = 0.0
        self.yoffset = 0.0
        self.fwhm = 0.0
        self.cax = 0.0
    
    def get_positions(self):
        self.position = (self.original_positions + self.xoffset)*self.xscale           
        return self.position
    def get_profilevalues(self):
        self.profilevalues = (self.original_values + self.yoffset)*self.yscale 
        return self.profilevalues
    
    def get_valfromposition(self,x): 
        if(self.position[1]-self.position[0])>0:
            cs = CubicSpline(self.get_positions(), self.get_profilevalues())
        else:
            cs = CubicSpline(np.flipud(self.get_positions()),np.flipud(self.get_profilevalues()))
        return cs(x)
    
    def get_cax(self, valref):
        _Y = self.get_profilevalues()
        _X = self.get_positions()
        _half_max = valref / 2.
        #find when function crosses line half_max (when sign of diff flips)
        #take the 'derivative' of signum(half_max - Y[])
        _d = np.sign(_half_max - np.array(_Y[0:-1])) - np.sign(_half_max - np.array(_Y[1:]))
        #plot(X[0:len(d)],d) #if you are interested
        #find the left and right most indexes
        _left_idx = np.where(_d > 0)[0]
        _right_idx = np.where(_d < 0)[-1]
        _arrx = _Y[_left_idx[0]:_left_idx[0]+2]
        _arry = _X[_left_idx[0]:_left_idx[0]+2]
        _f = interp1d(_arrx, _arry)
        xleft = _f(_half_max)
        _arrx = _Y[_right_idx[0]:_right_idx[0]+2]
        _arry = _X[_right_idx[0]:_right_idx[0]+2]
        _f = interp1d(_arrx, _arry)
        xright = _f(_half_max)
        self.fwhm = xright-xleft
        self.cax = xleft + (xright-xleft)/2
        return self.cax
    
    def get_maxval(self):
        return np.max(self.get_profilevalues())
    
    def get_FWHM(self):
        _cax = self.get_cax(np.max(self.get_profilevalues()))
        self.get_cax(self.get_valfromposition(_cax))
        return self.fwhm
    
    def get_L50(self):
        _cax = self.get_cax(np.max(self.get_profilevalues()))
        self.get_cax(self.get_valfromposition(_cax))
        return (self.cax - self.fwhm/2)
    
    def get_R50(self):
        _cax = self.get_cax(np.max(self.get_profilevalues()))
        self.get_cax(self.get_valfromposition(_cax))
        return (self.cax + self.fwhm/2)
    
    def center_cax(self):
        self.get_FWHM()
        self.xoffset = self.xoffset - self.cax/self.xscale
        
    def normalize_cax(self):
        self.get_FWHM()
        self.yscale = self.yscale/self.get_valfromposition(self.cax)
        
    def normalize_to_max(self):
        self.yscale = self.yscale/np.max(self.get_profilevalues())
    
    def profiles_difference(self, profile2):
        _minpos = np.max(np.array([np.min(self.get_positions()), np.min(profile2.get_positions())]))
        _maxpos = np.min(np.array([np.max(self.get_positions()), np.max(profile2.get_positions())]))
        _numpoints = 200
        _xarr = np.linspace(_minpos,_maxpos,_numpoints)
        _y1arr = np.array([self.get_valfromposition(x) for x in _xarr])
        _y2arr = np.array([profile2.get_valfromposition(x) for x in _xarr])
        _yarr = _y1arr-_y2arr
        _mat = np.array([_xarr,_yarr])
        diff_profile = Profile(_mat.transpose())
        return diff_profile
    
    def profiles_division(self, profile2):
        _minpos = np.max(np.array([np.min(self.get_positions()), np.min(profile2.get_positions())]))
        _maxpos = np.min(np.array([np.max(self.get_positions()), np.max(profile2.get_positions())]))
        _numpoints = 200
        _xarr = np.linspace(_minpos,_maxpos,_numpoints)
        _y1arr = np.array([self.get_valfromposition(x) for x in _xarr])
        _y2arr = np.array([profile2.get_valfromposition(x) for x in _xarr])
        _yarr = np.divide(_y1arr,_y2arr)
        _mat = np.array([_xarr,_yarr])
        diff_profile = Profile(_mat.transpose())
        return diff_profile
    
    