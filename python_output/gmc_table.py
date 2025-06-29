''' Output tables in a neat way'''
''' For paper, a few example lines, Latex format, decimal points aligned '''
''' For machine-readable version, full table, csv format,  decimal points aligned and fixed column widths'''
''' Fu-Heng Eric Liang 2024 Apr, ericfuhengliang@gmail.com '''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u

def get_digit_rough(value):
    sign_flag = 0
    if value < 0.:
        value = (-1) * value
        sign_flag = 1
    
    if value >= 1.0: digit = int(np.log10(value))+1
    else: digit = 1 # a zero is always there
    
    return digit # + sign_flag

def get_decimal_err(value):
    lowest_allowed = 2 # This is a bit arbitrary. Could change to 1 in principle.
    if round(value) >= lowest_allowed:
        return 0
    else:
        trial = np.format_float_positional(float("{0:.1g}".format(value)),trim='-') # to avoid potential scientific/exponential notation
        if '.' not in trial: # it must be less than lowest_allowed, eg 1
            return 1 # thus it needs one decimal digit, eg 1.3
        else:
            if (int(trial[-1]) < lowest_allowed):
                return len(trial.split('.')[1])+1
            else:
                return len(trial.split('.')[1])

def get_max_digit_both(array_raw, array_error_raw, select=None):
    if select is None: select = np.ones(len(array_raw),dtype=bool)
    max_value = -1
    max_error = -1
    max_sign = 0
    select_combine = select & np.isfinite(array_raw) & np.isfinite(array_error_raw)
    array = array_raw[select_combine]
    array_error = array_error_raw[select_combine]
    for i,v in enumerate(array):
        precision = get_decimal_err(array_error[i])
        value_rounded = float("{:.{}f}".format(v,precision))
        error_rounded = float("{:.{}f}".format(array_error[i],precision))
        value_digit = get_digit_rough(value_rounded)
        error_digit = get_digit_rough(error_rounded)
        if value_digit > max_value: max_value = value_digit
        if error_digit > max_error: max_error = error_digit
        if value_rounded < 0.: # return False from -0.0 < 0.
            max_sign = 1
    max_decimal_err = get_decimal_err(np.nanmin(array_error))
    # if np.min(array) < 0: max_sign = 1 # Mind the case of -0.0007 +- 0.1 --> 0.0 +- 0.1 without a sign
    # else: max_sign = 0
    return [max_value,max_error,max_decimal_err,max_sign]

def phrase(array, array_error, i, max_digit_both, flag=True, error_text=True):

    max_digit = max_digit_both[0]
    max_digit_err = max_digit_both[1]
    max_decimal_err = max_digit_both[2]
    max_sign = max_digit_both[3]
    max_decimal = max_decimal_err
    max_point_err = np.clip(max_decimal_err,0,1)
    max_point = max_point_err

    this_value = array[i]
    this_value_err = array_error[i]

    bad_flag = False
    if (~np.isfinite(this_value)) or (~np.isfinite(this_value_err) or (not flag)):
        this_value = 0
        this_value_err = 8
        bad_flag = True

    this_decimal_err = get_decimal_err(this_value_err)
    this_decimal = this_decimal_err
    this_point_err = np.clip(this_decimal_err,0,1)
    this_point = this_point_err
    value_rounded = float("{:.{}f}".format(this_value,this_decimal))
    error_rounded = float("{:.{}f}".format(this_value_err,this_decimal))
    this_digit = get_digit_rough(value_rounded)
    this_digit_err = get_digit_rough(error_rounded)

    if value_rounded == (-1)*value_rounded: this_value = np.absolute(this_value) # for the case of "-0.0", to get rid of "-"
    if this_value < 0: this_sign = 1
    else: this_sign = 0

    phrase_full = "{:>{front}.{pre}f}{xx:<{end}}".format(this_value, front=(max_sign+max_digit+this_point+this_decimal), pre=this_decimal, xx='',end=max_point+max_decimal-this_point-this_decimal)

    phrase_short = '$'+'\\phantom{'+ '-'*(max_sign-this_sign)+'0'*(max_digit-this_digit) +'}$'+'${:.{}f}'.format(this_value, this_decimal) + '\\phantom{'+'.'*(max_point-this_point) + '0'*(max_decimal-this_decimal)+'}'
    # if (this_sign == 0) and (max_sign == 1):
    #     phrase_short = '$\\phantom{-}$' + phrase_short

    if error_text:
        phrase_full = phrase_full+'+-'+"{:>{front}.{pre}f}{xx:<{end}}".format(this_value_err, front=max_digit_err+this_point_err+this_decimal_err, pre=this_decimal_err, xx='',end=max_point_err+max_decimal_err-this_point_err-this_decimal_err)
        
        phrase_short = phrase_short + ' \\pm ' + '\\phantom{'+ '0'*(max_digit_err-this_digit_err) +'}' + '{:.{}f}'.format(this_value_err, this_decimal_err) + '\\phantom{'+'.'*(max_point_err-this_point_err) + '0'*(max_decimal_err-this_decimal_err)+'}'

    if bad_flag:
        phrase_full = '-'.center(len(phrase_full))
        # print(repr(phrase_short))
        # print(phrase_short)
        if error_text: phrase_short = '-'.center(len(repr(phrase_short))-6)
        else: phrase_short = '-'.center(len(repr(phrase_short))-3)
    else:
        phrase_short = phrase_short+'$'
        pass

    return phrase_short+' & ', phrase_full+', '


col_number = 12
short_entry = 10 # number of (upper) entries in the table

# indfac = 4.0 # sqrt(pix_per_beam) = 4.0; applied in cprops bootstrap step to account for spatial correlation. Make sense or not?

inprop='/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits'
kinematic_table = '/Users/ericliang/n1387/work_pub/measurements/NGC1387_angmom_comparison.fits'
ingalpar='/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat'
resol_file = '/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt'
unres_file = '/Users/ericliang/n1387/work_pub/measurements/NGC1387-unresolved_bool.txt'
lines = open(ingalpar).readlines()

out_short = '/Users/ericliang/n1387/work_pub/plot/table_short.tex'
out_full = '/Users/ericliang/n1387/work_pub/plot/table_full.csv'

''' remove "non-clouds" '''
resol_ind = np.loadtxt(resol_file)
unres_ind = np.loadtxt(unres_file)
select_all = (resol_ind + unres_ind) > 0.5


catalog = fits.getdata(inprop)[select_all]
kine = fits.getdata(kinematic_table)[0]

# in unit of arcsec  # pixel size 0.04", resolution 0.15"
delta_ra = (catalog['mom1x'] - catalog['gal_cen_x']) * (-1) * round(catalog['degperpix'][0]*3600,3) 
delta_ra_error = np.ones(len(delta_ra)) * 0.04
delta_dec = (catalog['mom1y'] - catalog['gal_cen_y']) * round(catalog['degperpix'][0]*3600,3) 
delta_dec_error = np.ones(len(delta_dec)) * 0.04

omega_o = kine['OMEGA_O'][select_all]
omega_o_err = kine['ERR_OMEGA_O'][select_all]
theta_o = kine['THETA_O'][select_all]
theta_o_err = kine['ERR_THETA_O'][select_all]

v_sys = float(lines[16].strip())
delta_v = catalog['vpos'] - v_sys
delta_v_error = np.ones(len(delta_v)) * 0.4 # because spectral resolution ~4 km/s, going down by 1 order of mag

r_error = catalog['radrms_extrap_deconv_uc'] * catalog['radrms_extrap_deconv']# / indfac 
vrms_error = catalog['vrms_extrap_deconv'] * catalog['vrms_extrap_deconv_uc']# / indfac 

luminosity = catalog['lum_extrap'] / 1e4 # in unit of 1e4 K km/s pc^2
luminosity_err = luminosity * catalog['lum_extrap_uc']# / indfac

alpha_co = float(lines[18].strip())
mass = luminosity * alpha_co / 10 # in unit of 1e5 M_sun
mass_err = luminosity_err * alpha_co / 10

tmax_err = np.ones(len(mass_err)) * 1.1 # RMS=1.1K in flat-noise cube

r_gal_err = np.ones(len(mass_err)) * 4 # 0.04" = 4 pc

id_all = np.arange(len(mass_err),dtype=int)
short_flag_all = (id_all<short_entry) | (id_all==(len(id_all)-1))

max_all = []
max_all.append(get_max_digit_both(delta_ra, delta_ra_error))
max_all.append(get_max_digit_both(delta_dec, delta_dec_error))
# max_all.append(get_max_digit_both(delta_v, delta_v_error))
max_all.append(get_max_digit_both(catalog['vpos'], delta_v_error))
max_all.append(get_max_digit_both(catalog['radrms_extrap_deconv'], r_error, catalog['RESOLVE_SPATIAL'] > 0.5))
max_all.append(get_max_digit_both(catalog['vrms_extrap_deconv'], vrms_error,catalog['RESOLVE_SPECTRAL'] > 0.5))
max_all.append(get_max_digit_both(luminosity, luminosity_err))
max_all.append(get_max_digit_both(mass, mass_err))
max_all.append(get_max_digit_both(catalog['maxval'], tmax_err))
max_all.append(get_max_digit_both(omega_o, omega_o_err,catalog['RESOLVE_SPATIAL']*catalog['RESOLVE_SPECTRAL']> 0.5))
max_all.append(get_max_digit_both(theta_o, theta_o_err,catalog['RESOLVE_SPATIAL']*catalog['RESOLVE_SPECTRAL']> 0.5))
max_all.append(get_max_digit_both(catalog['r_gal'], r_gal_err))
max_all = np.array(max_all)

max_short = []
max_short.append(get_max_digit_both(delta_ra, delta_ra_error,short_flag_all))
max_short.append(get_max_digit_both(delta_dec, delta_dec_error,short_flag_all))
# max_short.append(get_max_digit_both(delta_v, delta_v_error,short_flag_all))
max_short.append(get_max_digit_both(catalog['vpos'], delta_v_error,short_flag_all))
max_short.append(get_max_digit_both(catalog['radrms_extrap_deconv'], r_error, (catalog['RESOLVE_SPATIAL']>0.5)&short_flag_all))
max_short.append(get_max_digit_both(catalog['vrms_extrap_deconv'], vrms_error,(catalog['RESOLVE_SPECTRAL']>.5)&short_flag_all))
max_short.append(get_max_digit_both(luminosity, luminosity_err,short_flag_all))
max_short.append(get_max_digit_both(mass, mass_err,short_flag_all))
max_short.append(get_max_digit_both(catalog['maxval'], tmax_err,short_flag_all))
max_short.append(get_max_digit_both(omega_o, omega_o_err,(catalog['RESOLVE_SPATIAL']*catalog['RESOLVE_SPECTRAL']>0.5)&short_flag_all))
max_short.append(get_max_digit_both(theta_o, theta_o_err,(catalog['RESOLVE_SPATIAL']*catalog['RESOLVE_SPECTRAL']>0.5)&short_flag_all))
max_short.append(get_max_digit_both(catalog['r_gal'], r_gal_err,short_flag_all))
max_short = np.array(max_short)


with open(out_short, 'w') as file_short:
    with open(out_full, 'w') as file_full:

        file_short.write('    \\begin{tabular}{'+'c'*col_number +'}\n')
        file_short.write('        \\hline\n')
        file_short.write('        \\hline\n')

        header_short = '        ID & RA & Dec. & $V_\\mathrm{bc}$ & $R_\\mathrm{c}$ & $\\sigma_\\mathrm{obs,los}$ & $L_\\mathrm{CO(2-1)}$ & $M_\\mathrm{gas}$ & $T_\\mathrm{max}$ & $\\omega_\\mathrm{obs}$ & $\\phi_\\mathrm{rot}$ & $R_\\mathrm{gal}$\\\\\n'
        header_full = '# ID, RA, Dec., V_bc, R_c, sigma_obs,los, L_CO(2-1), M_gas, T_max, omega_obs, phi_rot, R_gal\n'

        unit_short = '         & (h:m:s) & ($\\degr:\\arcmin:\\arcsec$) & (km~s$^{-1}$) & (pc) & (km~s$^{-1}$) & ($10^4$~K~km~s$^{-1}$~pc$^2$) & ($10^5$~M$_\odot$) & (K) & (km~s$^{-1}$~pc$^{-1}$) & (degree) & (pc)\\\\\n        \\hline\n'
        unit_full = '# , (h:m:s), (degree:arcmin:arcsec), (km s^-1), (pc), (km s^-1), (10^4 K km s^-1 pc^2), (10^5 M_\odot), (K), (km s^-1 pc^-1), (degree), (pc)\n'

        file_short.write(header_short + unit_short)
        file_full.write(header_full + unit_full)

        for i in range(len(omega_o)):

            catalog_this = catalog[i]
            short_flag = (i<short_entry) or (i==len(omega_o)-1)

            file_full.write(str(i+1).rjust(len(str(len(omega_o)))) + ', ')
            if short_flag: file_short.write( '        $\\phantom{'+'0'*(len(str(len(omega_o)))-len(str(i+1)))+'}'+ str(i+1)+'$ & ')

            # phrase_short,phrase_full = phrase(delta_ra,delta_ra_error,i,max_all[0],error_text=False)
            # file_full.write(phrase_full);
            # if short_flag: 
            #     phrase_short,phrase_full = phrase(delta_ra,delta_ra_error,i,max_short[0],error_text=False)
            #     file_short.write(phrase_short)

            # phrase_short,phrase_full = phrase(delta_dec,delta_dec_error,i,max_all[1],error_text=False)
            # file_full.write(phrase_full);
            # if short_flag: 
            #     phrase_short,phrase_full = phrase(delta_dec,delta_dec_error,i,max_short[1],error_text=False)
            #     file_short.write(phrase_short)

            ra_this = Angle(catalog_this['xpos'] * u.deg).hms
            phrase_tex = '$%d$:$%d$:$%.3f$' %(ra_this[0],ra_this[1],ra_this[2]) 
            phrase_csv = '%d:%d:%.3f' %(ra_this[0],ra_this[1],ra_this[2]) 
            file_full.write(phrase_csv+', ')
            if short_flag: 
                file_short.write(phrase_tex+' & ')

            dec_this = Angle(catalog_this['ypos'] * u.deg).signed_dms
            if dec_this[0] < 0: dec_sign = '-'
            else: dec_sign = '+'
            phrase_tex = '$'+dec_sign+'%d$:$%d$:$%.2f$' %(dec_this[1],dec_this[2],dec_this[3]) 
            phrase_csv = dec_sign+'%d:%d:%.2f' %(dec_this[1],dec_this[2],dec_this[3]) 
            file_full.write(phrase_csv+', ')
            if short_flag: 
                file_short.write(phrase_tex+' & ')

            # phrase_short,phrase_full = phrase(delta_v,delta_v_error,i,max_all[2],error_text=False)
            # file_full.write(phrase_full);
            # if short_flag: 
            #     phrase_short,phrase_full = phrase(delta_v,delta_v_error,i,max_short[2],error_text=False)
            #     file_short.write(phrase_short)

            phrase_short,phrase_full = phrase(catalog['vpos'], delta_v_error,i,max_all[2],error_text=False)
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full = phrase(catalog['vpos'], delta_v_error,i,max_short[2],error_text=False)
                file_short.write(phrase_short)


            phrase_short,phrase_full = phrase(catalog['radrms_extrap_deconv'],r_error,i,max_all[3],catalog_this['RESOLVE_SPATIAL']>0.5)
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(catalog['radrms_extrap_deconv'],r_error,i,max_short[3],catalog_this['RESOLVE_SPATIAL']>0.5)
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(catalog['vrms_extrap_deconv'],vrms_error,i,max_all[4],catalog_this['RESOLVE_SPECTRAL']>0.5)
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(catalog['vrms_extrap_deconv'],vrms_error,i,max_short[4],catalog_this['RESOLVE_SPECTRAL']>0.5)
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(luminosity,luminosity_err,i,max_all[5])
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(luminosity,luminosity_err,i,max_short[5])
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(mass,mass_err,i,max_all[6])
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(mass,mass_err,i,max_short[6])
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(catalog['maxval'],tmax_err,i,max_all[7],error_text=False)
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(catalog['maxval'],tmax_err,i,max_short[7],error_text=False)
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(omega_o,omega_o_err,i,max_all[8],catalog_this['RESOLVE_SPATIAL']*catalog_this['RESOLVE_SPECTRAL']>0.5)
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(omega_o,omega_o_err,i,max_short[8],catalog_this['RESOLVE_SPATIAL']*catalog_this['RESOLVE_SPECTRAL']>0.5)
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(theta_o,theta_o_err,i,max_all[9],catalog_this['RESOLVE_SPATIAL']*catalog_this['RESOLVE_SPECTRAL']>0.5)
            file_full.write(phrase_full);
            if short_flag: 
                phrase_short,phrase_full=phrase(theta_o,theta_o_err,i,max_short[9],catalog_this['RESOLVE_SPATIAL']*catalog_this['RESOLVE_SPECTRAL']>0.5)
                file_short.write(phrase_short)

            phrase_short,phrase_full=phrase(catalog['r_gal'], r_gal_err,i,max_all[10],error_text=False)
            file_full.write(phrase_full[:-2]);
            if short_flag: 
                phrase_short,phrase_full=phrase(catalog['r_gal'], r_gal_err,i,max_short[10],error_text=False)
                file_short.write(phrase_short[:-2])


            file_full.write('\n');
            if short_flag: file_short.write('\\\\\n')

            if i == short_entry:
                file_short.write(('        '+'... & '*col_number)[:-3]+'\\\\\n')

        file_short.write('        \\hline\n')
        file_short.write('        \\hline\n')
        file_short.write('    \\end{tabular}')




# def get_max_digit(array):
#     return np.max([get_digit_rough(np.nanmax(array)),get_digit_rough(np.nanmin(array))])

            # max_digit = get_max_digit(delta_ra)
            # max_point = 1; max_decimal = 2 # because pixel size is 0.04", resolution 0.15"
            # this_digit = get_digit_rough(delta_ra[i])
            # this_point = 1; this_decimal = 2

            # file_full.write('{0: >{1}.{2}f}'.format(delta_ra[i], (max_digit+max_point+max_decimal), this_decimal) +' , ')
            # if short_flag: file_short.write('$\\phantom{'+ '0'*(max_digit+max_point+max_decimal-this_digit-this_point-this_decimal) +'}' + '{0:.{1}f}'.format(delta_ra[i], this_decimal) +'$ & ')

            # max_digit = get_max_digit(delta_dec)
            # this_digit = get_digit_rough(delta_dec[i])

            # file_full.write('{0: >{1}.{2}f}'.format(delta_dec[i], (max_digit+max_point+max_decimal), this_decimal) +' , ')
            # if short_flag: file_short.write('$\\phantom{'+ '0'*(max_digit+max_point+max_decimal-this_digit-this_point-this_decimal) +'}' + '{0:.{1}f}'.format(delta_dec[i], this_decimal) +'$ & ')

            # max_digit = get_max_digit(delta_v)
            # max_point = 1; max_decimal = 1 # because spectral resolution ~4 km/s, going down by 1 order of mag
            # this_digit = get_digit_rough(delta_v[i])
            # this_point = 1; this_decimal = 1

            # file_full.write('{0: >{1}.{2}f}'.format(delta_v[i], (max_digit+max_point+max_decimal), this_decimal) +' , ')
            # if short_flag: file_short.write('$\\phantom{'+ '0'*(max_digit+max_point+max_decimal-this_digit-this_point-this_decimal) +'}' + '{0:.{1}f}'.format(delta_v[i], this_decimal) +'$ & ')

