/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Rui Hu, while the main part of GADGET4 N-body/SPH code is 
 * \copyright	developed by Volker Springel. Copyright (C) 2014-2020 by 
 * \copyright	Volker Springel (vspringel@mpa-garching.mpg.de) and all
 * \copyright   contributing authors.
 *******************************************************************************/

/*! \neutrino.cc 
* 
*	\brief module for neutrino evolution
*/

#ifdef NEUTRINO

#include <gsl/gsl_rng.h>
#include <stdio.h>

#include "../neutrino/neutrino.h"
#include "../data/allvars.h"
#include "../data/mymalloc.h"
#include "../data/dtypes.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/** \brief Compute the evolution of neutrino.
* 
*/

void frsfr::InitNu(void) 
{ 
	All.H0 = All.Hubble * All.HubbleParam * 1000. / (1E6 * 3.0857E16);
	All.Rhocr = (cc * 3.0 * All.H0 * All.H0) / (8.0 * M_PI * Gr);
    All.unittrans = pow(kb, 4) / ((pow(hbar, 3)) * (pow(c, 3)) * 2 * M_PI * M_PI);
    Check_Neutrinos();

}

void frsfr::Check_Neutrinos(void) 
{
    switch(All.LeptonAsymmetry)
    {
        case 2:
           All.Xi_2 = 0;
           All.Xi_1 = 0;
           break;

        case 1: 
           All.Xi_2 = Cal_Xi2(All.Xi_3);
           All.Xi_1 = Cal_Xi1(All.Xi_2, All.Xi_3);
           break;        

        case 0: 
           All.Xi2 = All.Xi_3;
           All.Xi1 = All.Xi_3;
           break;

        default:
           break;
    }

    switch(All.MassHierarchy)
    {
        case 0: /* Normal */
          All.Mass_2 = sqrt(All.Mass_1 * All.Mass_1 + 7.59e-5);
          All.Mass_3 = sqrt(All.Mass_2 * All.Mass_2 + 2.32e-3);
          break;

        case 1: /* Inverted */
          All.Mass_2 = sqrt(All.Mass_1 * All.Mass_1 + 2.32e-3);
          All.Mass_3 = sqrt(All.Mass_2 * All.Mass_2 + 7.59e-5);
          break;

        case 2: /* Identical */
          All.Mass_2 = All.Mass_1;
          All.Mass_3 = All.Mass_1;
          break;

        case 3: /* for only one sterile/active neutrino */
          All.Mass_2 = 0;
          All.Mass_3 = All.Mass_1;
          All.Mass_1 = 0;
          break;

        default:
          break;
    }
    
    printf("Mass1: %f Mass2: %f Mass3: %f Xi_1: %f Xi_2: %f Xi_3: %f \n", All.Mass_1, All.Mass_2, All.Mass_3, All.Xi_1, All.Xi_2,
           All.Xi_3);

    switch(All.ExpanOn)
    {
        case 2: /* for only one sterile/active neutrino */
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, All.Mass_3, All.Xi_3);
          All.Omega_Nu0_Frstr     = neutrino_integration(1.0, 0., 0.);
          break;

        case 1:
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, All.Mass_3, All.Xi_3) + neutrino_integration(1.0, All.Mass_2, All.Xi_2) +
                                    neutrino_integration(1.0, All.Mass_1, All.Xi_1);
          All.
          break;

        case 0:
          All.Omega_Nu0_Expansion = neutrino_integration(1.0, 0., 0.) * 3;
          break;

        default:
          break;
    }

}
