# -*- coding: utf-8 -*-

"""Constants for unit tests."""

import os

from prodec.utils import std_amino_acids

TESTS_FLDR = os.path.dirname(__file__)
SRC_FLDR = os.path.join(os.path.dirname(TESTS_FLDR), 'src', 'prodec')


DFLT_DATA = os.path.join(SRC_FLDR, 'data.json')

DFLT_AA = list('ACDEFGHIKLMNPQRSTVWY')

CSTM_DATA_NULL = """{
      "ID": "CSTM_TESTS", "Type": "Linear", "Name": "TEST",
      "Authors": null, "Year": null, "Journal": null,
      "DOI": null, "PMID": null, "Patent": null,
      "Scales": { "Names": [ "v1", "v2", "v3" ],
         "Values": { A_value, C_value, 
                    D_value, E_value, F_value, 
                    G_value, H_value, I_value, 
                    K_value, L_value, M_value, 
                    N_value, P_value, Q_value, 
                    R_value, S_value, T_value, 
                    V_value, W_value, Y_value
         }
      }
   }"""

DFLT_SEQ = ''.join(std_amino_acids)
STUPID_SEQ = 'AZERTYUIOP'
