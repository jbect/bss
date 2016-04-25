/******************************************************************************
 *                                                                            *
 * Copyright Notice                                                           *
 *                                                                            *
 *    Copyright (C) 2016 CentraleSupelec                                      *
 *                                                                            *
 *    Author:  Julien Bect  <julien.bect@centralesupelec.fr>                  *
 *                                                                            *
 *    Address questions, bug reports, feature requests, or any other          *
 *    correspondance related to BSS to kriging-help@lists.sourceforge.net     *
 *                                                                            *
 * Copying Permission Statement                                               *
 *                                                                            *
 *    This file is part of BSS (https://github.com/jbect/bss).                *
 *                                                                            *
 *    BSS is free software; you can redistribute it and/or modify it under    *
 *    the terms of the  GNU Lesser General Public License  as published by    *
 *    the Free Software Foundation;  either version 2.1 of the License, or    *
 *    (at your option) any later version.                                     *
 *                                                                            *
 *    BSS is distributed  in the hope that it will  be useful, but WITHOUT    *
 *    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY    *
 *    or FITNESS  FOR A  PARTICULAR PURPOSE.  See the  GNU  Lesser General    *
 *    Public License for more details.                                        *
 *                                                                            *
 *    You should  have received  a copy  of the  GNU Lesser General Public    *
 *    License along with BSS;  if not, see <http://www.gnu.org/licenses/>.    *
 *                                                                            *
 ******************************************************************************/

#include "mex.h"

#define X_IN  prhs[0]
#define COUNT prhs[1]
#define X_OUT plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t m_in, d, m1, m2, k, i_in, i_out, j, nb_copies, m_out;
    double *x_in, *x_out, *c, *q_in, *q_out, s;
    
    if(nrhs != 2)
        mexErrMsgTxt("Incorrect number of input arguments.");
    
    if(nlhs > 1)
        mexErrMsgTxt("Incorrect number of output arguments.");
    
    if(mxGetNumberOfDimensions(X_IN) != 2)
        mexErrMsgTxt("x_in must be a matrix.");
    
    m_in = mxGetM(X_IN);
    d = mxGetN(X_IN);
    
    if(mxGetNumberOfDimensions(COUNT) != 2)
        mexErrMsgTxt("c must be a vector.");
    
    m1 = mxGetM(COUNT);
    m2 = mxGetN(COUNT);
    
    if(!(((m1 == m_in) && (m2 == 1)) || ((m1 == 1) && (m2 == m_in))))
        mexErrMsgTxt("c must be a vector, with length "
                "equal to the number of rows in x_in");
    
    c = mxGetPr(COUNT);
    
    s = 0;
    for(i_in = 0; i_in < m_in; i_in++)
        s += c[i_in];
    m_out = (size_t) s;
    
    X_OUT = mxCreateDoubleMatrix(m_out, d, mxREAL);
    
    x_in  = mxGetPr(X_IN);
    x_out = mxGetPr(X_OUT);
    
    for(i_in = 0, i_out = 0; i_in < m_in; i_in++, i_out += nb_copies)
    {
        nb_copies = (size_t) c[i_in];
        
        if(nb_copies > 0)
        {
            q_in  = x_in  + i_in;
            q_out = x_out + i_out;
            
            /* loop over dimensions */
            for(k = 0; k < d; k++, q_in += m_in, q_out += m_out)
                for(j = 0; j < nb_copies; j++)
                    q_out[j] = q_in[0];
        }
    }
}