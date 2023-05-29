
def find_J(a_j, ci, gammastar):
    """
    Function finds J (light limited electron flux)
    Equation comes from Morfopolous PhD Chapter 2, Equation 2.16 for Aj

    Args: # get the input units right as well as the output
        a_j (float): light limmited gross assimilation (umol m-2 s-1)
        ci (float): internal CO2 concentration (Pa)
        gammastar (float): compensation point in the absance of dark resperation (Pa)
    
    Returns:
        J (float): light limited electron flux (electron flux used for C and O2 assimilation) (umols m-2 s-1)
    """
    J = 4 * a_j * ((ci + 2*gammastar) / (ci - gammastar))
    return J

def find_Jv(vcmax, ci, gammastar, kmm):
    """
    Function finds Jv (electron flux for Rubisco limmited assimilation)
    Formula comes from Morfopolous PhD Chapter 2, Equation 2.26*

    Args: #need to get the units
        vcmax (float): maximum carboxylation capacity (umols m-2 s-1)
        ci (float): light limmited gross assimilation (ubar)
        gammastar (float): compensation point in the absance of dark resperation (ubar)
        kmm (float): Michalis Menten coefficient for C, when there is competition for Rubisco by O2 and C (ubar)
    
    Returns:
        Jv (float): Rubisco limited electron flux (umols m-2 s-1)
    """    
    # kmm = kc_dashed 
    #Jv = 4 * vcmax * ((ci - gammastar)/(ci + kc_dashed))

    #input comes in in Pa , assume Vcmax is correct for now
    # convert Pa to ubar
    ci_ubar = ci *10
    gammastar_ubar = gammastar *10
    kmm_ubar = kmm *10


    Jv = 4 * vcmax * ((ci_ubar + 2*gammastar_ubar)/(ci_ubar + kmm_ubar)) #from the new physiologist
    return Jv