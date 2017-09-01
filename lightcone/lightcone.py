from __future__ import absolute_import, division, print_function
import numpy as np
import warnings
# from .due import due, Doi

__all__ = ["tao_impl_angle_beta"]


import astropy
import astropy.units as u

def convert_redshift_to_comoving_distance(redshifts,
                                          cosmo=None,
                                          distance_units=None):
    r"""
    Returns co-moving distances for a list of redshifts

    Parameters:
    -----------
    
    redshifts: double, array or scalar
               List of redshifts
    cosmo    : An astropy cosmology object. default is Planck 2015
               Sets up the cosmology for calculating the co-moving distances
    distance_units : astropy units object, units for the co-moving distance(s)
       
    Returns:
    --------
    com_dist : double, array or scalar. Same length as the input redshifts
               Returned in `distance_units` (if specified); otherwise, returned
               in Mpc/h units. 
    """

    default_units = u.Mpc
    
    if cosmo is None:
        from astropy.cosmology import Planck15
        cosmo = Planck15
        msg = "cosmology is not set. Using Planck 2015 cosmology = {0}"\
              .format(cosmo)
        warnings.warn(msg)
    else:
        if not isinstance(cosmo, astropy.cosmology):
            msg = '{The cosmology object parameter must be an instance of '\
                  '`astropy.cosmology`}'
            raise ValueError(msg)

    if distance_units is not None:
        if not isinstance(distance_units, u.Unit):
            msg = 'distance units parameter = {0} must be an instance of '\
                  '`astropy.units`.'.format(distance_units)
            raise u.UnitsError(msg)
    else:
        distance_units = default_units


    # calculate all the distances
    distances = cosmo.comoving_distance(redshifts)
    distances.to(default_units)


    # calculate 1/h
    H0 = cosmo.H(0)
    default_hubble_units = (u.km/u.s)/u.Mpc
    hundredkm_per_s_per_Mpc = 100.0 * default_hubble_units
    little_h = cosmo.H(0).to(default_hubble_units)/hundredkm_per_s_per_Mpc
    print("H0 = {0} little h = {1}".format(H0, little_h))

    # convert to co-moving 1/h units
    distances = distances * little_h
    print("distances = {0}".format(distances))

    # Now return in the requested units
    # (set to Mpc/h by default)
    return distances.to(distance_units)


        
# # Use duecredit (duecredit.org) to provide a citation to relevant work to
# # be cited. This does nothing, unless the user has duecredit installed,
# # And calls this with duecredit (as in `python -m duecredit script.py`):
# due.cite(Doi("10.1167/13.9.30"),
#          description="Template project for small scientific Python projects",
#          tags=["reference-implementation"],
#          path='lightcone')

def tao_paper_solution_angle_beta(min_ra=10.0, max_ra=20.0,
                                  min_dec=30.0, max_dec=35.0,
                                  zmin=0.0, zmax=2.0,
                                  boxsize=[500.0, 500.0, 500.0]):
    r"""
    Returns the angle, :math:`\beta`, to construct an unique lightcone

    The routine here is an attempt at creating the algorithm presented
    in the TAO code paper (http://adsabs.harvard.edu/abs/2016ApJS..223....9B)

    Parameters:
    -----------

    min_ra  : double, defaut=10.0 degrees. Must be in range [0.0, 360.0]
              Minimum value of Right Ascension for the lightcone

    max_ra  : double, default=20.0 degrees. Must be in range [0.0, 360.0]
              Maximum value of Right Ascension for the lightcone

    min_dec : double, default=30.0 degrees, Must be in range [-90.0, 90.0]
              Minimum value of Declination for the lightcone

    min_dec : double, default=30.0 degrees, Must be in range [-90.0, 90.0]
              Maximum value of Declination for the lightcone

    zmin    : double, default=0.0. 
              Minimum redshift cut for the lightcone

    zmax    : double, default=2.0
              Maximum redshift cut for the lightcone

    boxsize : double, or array of 3 doubles, units=Mpc/h. default = 500 Mpc/h. 
              The periodic boxsize in each of the 3 dimensions.
    
    Returns:
    --------

    beta : double, units=degrees
           The angle by which the lightcone needs to start off such that an
           unique lightcone solution can be generated. 
            
           
    
    .. note : The solution here might be different from the one presented in 
              the TAO paper (http://adsabs.harvard.edu/abs/2016ApJS..223....9B)
    
    
    
    """
    



def tao_impl_angle_beta(min_ra=10.0, max_ra=20.0,
                        min_dec=30.0, max_dec=35.0,
                        zmin=0.0, zmax=2.0,
                        cosmo=None,
                        boxsize=500.0 * u.Mpc):
                        
    r"""
    Returns the angle, :math:`\beta`, to construct an unique lightcone

    The routine here is a direct translation from C++ in the TAO code-base
    to python.

    Parameters:
    -----------

    min_ra  : double, defaut=10.0 degrees. Must be in range [0.0, 360.0]
              Minimum value of Right Ascension for the lightcone

    max_ra  : double, default=20.0 degrees. Must be in range [0.0, 360.0]
              Maximum value of Right Ascension for the lightcone

    min_dec : double, default=30.0 degrees, Must be in range [-90.0, 90.0]
              Minimum value of Declination for the lightcone

    min_dec : double, default=30.0 degrees, Must be in range [-90.0, 90.0]
              Maximum value of Declination for the lightcone

    zmin    : double, default=0.0. 
              Minimum redshift cut for the lightcone

    zmax    : double, default=2.0
              Maximum redshift cut for the lightcone

    cosmo   : astropy cosmology object. default None

    boxsize : double, or array of 3 doubles, units=Mpc/h. default = 500 Mpc/h. 
              The periodic boxsize in each of the 3 dimensions.
    
    Returns:
    --------

    beta : double, units=degrees
           The angle by which the lightcone needs to start off such that an
           unique lightcone solution can be generated. 
            
           
    
    .. note : The solution here might be different from the one presented in 
              the TAO paper (http://adsabs.harvard.edu/abs/2016ApJS..223....9B)
    
    
    """
    max_redshift_allowed = 100.0
    
    # Input validation
    if min_ra < 0.0 or max_ra > 360:
        msg = 'Right Ascension (RA) must be between [0.0, 360.0]. The input '\
              'RA min, max values are = {0}, {1}'.format(min_ra, max_ra)
        raise ValueError(msg)

    if min_dec < -90.0 or max_dec > 90.0:
        msg = 'Declination (DEC) must be between [-90.0, 90.0]. The input '\
              'DEC min, max values are = {0}, {1}'.format(min_dec, max_dec)
        raise ValueError(msg)

    # Now attach degrees units
    print("Before switching to radians: min ra = {0}".format(min_ra))
    min_ra = (min_ra * u.deg).to(u.rad)
    max_ra = (max_ra * u.deg).to(u.rad)
    min_dec = (min_dec * u.deg).to(u.rad)
    max_dec = (max_dec * u.deg).to(u.rad)
    print("After switching to radians: min ra = {0}".format(min_ra))

    if zmin < 0.0 or zmax > max_redshift_allowed:
        msg = 'Redshift (z) must be between [0.0, {0}]. The input '\
              'z min, max values are = {1}, {2}'.format(max_redshift_allowed,
                                                        zmin,
                                                        zmax)
        raise ValueError(msg)

    units = None
    if isinstance(boxsize, u.Quantity):
        units = boxsize.unit
    
    d1, d0 = convert_redshift_to_comoving_distance([zmax, zmin],
                                                   cosmo=cosmo,
                                                   distance_units=units)
    # If boxsize did not have units, convert to
    # the units of the co-moving distance
    if units is None:
        boxsize = boxsize * d1.unit

    ra_diff = max_ra - min_ra
    if (d1 - d0 * np.cos(ra_diff)) <= boxsize:
        # all angles are in radians -> convert to degrees
        return 0.0 - min_ra.to(deg)
    
    # Use Ridder's method to find the optimal angle for unique cones
    # auto res = hpc::ridders(
    #     [ra, d0, d1, b]( double x )
    #     {  
    #         double phi = ra + x;
    #         return b - d1*(cos( x ) - sin( x )/tan( phi ));
    # },
    #     0.5*M_PI,
    #     0.0
    # );
    # if( res != std::numeric_limits<double>::max() )
    # return res - lc.min_ra();
    # else
    # return boost::none;
    def _func(x, ra_diff, d0, d1, boxsize):
        x_in_rad = x * u.rad
        for name, a in zip(['x', 'ra diff', 'd0', 'd1', 'boxsize'],
                           [x_in_rad, ra_diff, d0, d1, boxsize]):
            print("{0} = {1}".format(name, a))
            
        phi = ra_diff + x_in_rad
        res =  boxsize - d1 * (np.cos(x) - np.sin(x)/np.tan(phi))
        print("res = {0}".format(res))
        return res
        
    method = 'ridder'
    if method == 'ridder':
        from scipy.optimize import ridder as solver
    else:
        from scipy.optimize import brentq as solver

    beta, root_obj = solver(_func, 0.0, 0.5*np.pi,
                            maxiter=100,
                            full_output=True,
                            disp=True,
                            args=(ra_diff, d0, d1, boxsize))
                  
    print("Root object = {0}".format(root_obj))
    print("Solved angle = {0}. Converged = {1}"
          .format(beta, root_obj.Converged))
