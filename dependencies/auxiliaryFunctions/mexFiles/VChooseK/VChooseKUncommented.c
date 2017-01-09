#include "mex.h"
#include <math.h>

#ifndef mwSize
#define mwSize  int
#define mwIndex int
#define MWSIZE_MAX 2147483647UL
#define MWSIZE_MAX 2147483647UL
#endif

#define WARN_LIMIT   500000000L 
#define ERROR_LIMIT 1000000000L  

mwSize GetOutputLen(mwSize n, mwSize k, double *C_double);
void BadInputTypeError(void);

void Elem2_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem2_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem2_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem2_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void Elem3_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem3_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem3_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem3_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void Elem4_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem4_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem4_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem4_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void Elem5_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY);
void Elem5_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY);
void Elem5_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY);
void Elem5_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY);

void ElemK_1Byte(INT8_T  *Xp, mwSize nX, INT8_T  *Yp, mwSize nY, mwSize K);
void ElemK_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY, mwSize K);
void ElemK_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY, mwSize K);
void ElemK_8Byte(double  *Xp, mwSize nX, double  *Yp, mwSize nY, mwSize K);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize nX, nY, K;
  const mxArray *X;
  double nY_double;
  int ElementSize;
  
  if (nrhs != 2) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput", "2 inputs required.");
  }
  if (nlhs > 1) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput", "1 output allowed.");
  }
  
  X = prhs[0];
  if (!mxIsNumeric(X) && !mxIsChar(X) & !mxIsLogical(X)) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput",
                       "Input array must be numerical, CHAR or LOGICAL.");
  }
  if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:BadNInput",
                       "Input K must be a numerical scalar.");
  }
  
  nX = mxGetNumberOfElements(X);
  K  = (mwSize) floor(mxGetScalar(prhs[1]) + 0.5);  
  
  if (nX == 0) {
     plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     return;
  }
  
  if (K <= 1) {
     if (K == 1) { 
        plhs[0] = mxDuplicateArray(X);
        mxSetM(plhs[0], nX);
        mxSetN(plhs[0], 1);
     } else if (K <= 0) {  
       plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     }
     return;
  }
  
  if (K >= nX) {
     if (K == nX) {
        plhs[0] = mxDuplicateArray(X);
        mxSetM(plhs[0], 1);
        mxSetN(plhs[0], nX);
     } else {
        mexWarnMsgIdAndTxt("JSimon:VChooseK:ShortInput",
                           "Too short input: K > N.");
        plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     }
     return;
  }
  
  nY          = GetOutputLen(nX, K, &nY_double);
  ElementSize = mxGetElementSize(X);
  nY_double  *= ElementSize;
  if (nY_double > WARN_LIMIT) {
     if (nY_double < ERROR_LIMIT && nY_double < MWSIZE_MAX) {
        mexWarnMsgIdAndTxt("JSimon:VChooseK:LargeOutput",
                           "Output will be large and take a long time.");
     } else {
        mexErrMsgIdAndTxt("JSimon:VChooseK:TooLargeOutput",
                          "Output would be too large.");
     }
  }
  
  plhs[0] = mxCreateNumericMatrix(nY, K, mxGetClassID(X), mxREAL);
  
  switch (K) {
    case 2: 
      switch (ElementSize) {
        case 8:  Elem2_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem2_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem2_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem2_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;
        
    case 3:  
      switch (ElementSize) {
        case 8:  Elem3_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem3_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem3_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem3_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;
        
    case 4:  
      switch (ElementSize) {
        case 8:  Elem4_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem4_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem4_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem4_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;
      
    case 5: 
      switch (ElementSize) {
        case 8:  Elem5_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY);
                 break;
        case 4:  Elem5_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY);
                 break;
        case 2:  Elem5_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY);
                 break;
        case 1:  Elem5_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY);
                 break;
        default: BadInputTypeError();
      }
      break;

    default:  
      switch (ElementSize) {
        case 8:  ElemK_8Byte(mxGetPr(X), nX, mxGetPr(plhs[0]), nY, K);
                 break;
        case 4:  ElemK_4Byte((INT32_T *) mxGetData(X), nX,
                             (INT32_T *) mxGetData(plhs[0]), nY, K);
                 break;
        case 2:  ElemK_2Byte((INT16_T *) mxGetData(X), nX,
                             (INT16_T *) mxGetData(plhs[0]), nY, K);
                 break;
        case 1:  ElemK_1Byte((INT8_T *) mxGetData(X), nX,
                             (INT8_T *) mxGetData(plhs[0]), nY, K);
                 break;
        default: BadInputTypeError();
      }
  }
  
  return;
}

void BadInputTypeError(void)
{
  mexErrMsgIdAndTxt("JSimon:VChooseK:BadInputType",
                    "Input must have 1, 2, 4 or 8 bytes per element.");
}

mwSize GetOutputLen(mwSize n, mwSize k, double *C_double)
{
  mwSize i, ai, bi;
  double ad, bd;
  
  bi = n - k;
  ai = bi + 1;
  ad = (double) ai;
  bd = (double) bi;
  for (i = 2; i <= k; i++) {
     ai += (ai * bi) / i;
     ad += (ad * bd) / i;
  }
  
  *C_double = ad;
  
  return (ai);
}

void Elem2_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  double *Xq, *Xf, *Y1, *Y2;
  
  Y1 = Yp;             
  Y2 = Y1 + nY;     
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {  
      *Y1++ = *Xp;     
      *Y2++ = *Xq;     
    }
  }

  return;
}

void Elem2_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  INT32_T *Xq, *Xf, *Y1, *Y2;

  Y1 = Yp;             
  Y2 = Y1 + nY;        
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {
      *Y1++ = *Xp;     
      *Y2++ = *Xq;     
    }
  }

  return;
}

void Elem2_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  INT16_T *Xq, *Xf, *Y1, *Y2;

  Y1 = Yp;
  Y2 = Y1 + nY;        
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {
      *Y1++ = *Xp;     
      *Y2++ = *Xq;     
    }
  }

  return;
}

void Elem2_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  
  INT8_T *Xq, *Xf, *Y1, *Y2;

  Y1 = Yp; 
  Y2 = Y1 + nY;        
  for (Xf = Xp + nX - 1; Xp < Xf; Xp++) {
    for (Xq = Xp + 1; Xq <= Xf; Xq++) {
      *Y1++ = *Xp;     
      *Y2++ = *Xq;     
    }
  }

  return;
}

void Elem3_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  double *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           
  Y2  = Y1 + nY;   
  Y3  = Y2 + nY;   
  X3f = Xp + nX;   
  X2f = X3f - 1;    
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {         
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  
        *Y2++ = *Xq;
        *Y3++ = *Xr;  
      }
    }
  }

  return;
}

void Elem3_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  INT32_T *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           
  Y2  = Y1 + nY;   
  Y3  = Y2 + nY;   
  X3f = Xp + nX;   
  X2f = X3f - 1;     
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {          
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  
        *Y2++ = *Xq;
        *Y3++ = *Xr;  
      }
    }
  }
  
  return;
}


void Elem3_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  INT16_T *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;           
  Y2  = Y1 + nY;   
  Y3  = Y2 + nY;   
  X3f = Xp + nX;   
  X2f = X3f - 1;     
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {    
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  
        *Y2++ = *Xq;
        *Y3++ = *Xr;  
      }
    }
  }
  
  return;
}


void Elem3_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  INT8_T *Xq, *Xr, *X2f, *X3f, *Y1, *Y2, *Y3;
  
  Y1  = Yp;      
  Y2  = Y1 + nY;     
  Y3  = Y2 + nY;     
  X3f = Xp + nX;    
  X2f = X3f - 1;      
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {     
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        *Y1++ = *Xp;  
        *Y2++ = *Xq;
        *Y3++ = *Xr;  
      }
    }
  }
  
  return;
}

void Elem4_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  double *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;              
  Y2  = Y1 + nY;      
  Y3  = Y2 + nY;      
  Y4  = Y3 + nY;      
  X4f = Xp + nX;      
  X2f = X4f - 2;        
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {      
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp; 
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs; 
        }
      }
    }
  }

  return;
}

void Elem4_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  INT32_T *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;            
  Y2  = Y1 + nY;         
  Y3  = Y2 + nY;        
  Y4  = Y3 + nY;        
  X4f = Xp + nX;
  X2f = X4f - 2;        
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {       
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;  
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs;  
        }
      }
    }
  }
  
  return;
}

void Elem4_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  INT16_T *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;            
  Y2  = Y1 + nY;      
  Y3  = Y2 + nY;      
  Y4  = Y3 + nY;      
  X4f = Xp + nX;      
  X2f = X4f - 2;        
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {      
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs; 
        }
      }
    }
  }
  
  return;
}

void Elem4_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  INT8_T *Xq, *Xr, *Xs, *X2f, *X4f, *Y1, *Y2, *Y3, *Y4;
  
  Y1  = Yp;              
  Y2  = Y1 + nY;      
  Y3  = Y2 + nY;      
  Y4  = Y3 + nY;      
  X4f = Xp + nX;      
  X2f = X4f - 2;        
  
  for ( ; Xp <= X2f; Xp++) {
    for (Xq = Xp + 1; Xq < X2f; Xq++) {
      for (Xr = Xq + 1; Xr <= X4f; Xr++) {    
        for (Xs = Xr + 1; Xs < X4f; Xs++) {
          *Y1++ = *Xp;  
          *Y2++ = *Xq;
          *Y3++ = *Xr;
          *Y4++ = *Xs;  
        }
      }
    }
  }
  
  return;
}

void Elem5_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY)
{
  double *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
         *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               
  Y2  = Y1 + nY;       
  Y3  = Y2 + nY;       
  Y4  = Y3 + nY;       
  Y5  = Y4 + nY;      
  X5f = Xp + nX;      
                          
  X3f = X5f - 2;    
                          
  X1f = X3f - 2;    
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp; 
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt; 
          }
        }
      }
    }
  }

  return;
}

void Elem5_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY)
{
  INT32_T *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
          *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               
  Y2  = Y1 + nY;       
  Y3  = Y2 + nY;       
  Y4  = Y3 + nY;       
  Y5  = Y4 + nY;       
  X5f = Xp + nX;       
  X3f = X5f - 2;         
  X1f = X3f - 2;         
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp; 
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  
          }
        }
      }
    }
  }

  return;
}

void Elem5_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY)
{
  INT16_T *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
          *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               
  Y2  = Y1 + nY;       
  Y3  = Y2 + nY;       
  Y4  = Y3 + nY;       
  Y5  = Y4 + nY;       
  X5f = Xp + nX;       
  X3f = X5f - 2;         
  X1f = X3f - 2;        
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp;
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  
          }
        }
      }
    }
  }
  
  return;
}

void Elem5_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY)
{
  INT8_T *Xq, *Xr, *Xs, *Xt, *X5f, *X3f, *X1f,
         *Y1, *Y2, *Y3, *Y4, *Y5;
  
  Y1  = Yp;               
  Y2  = Y1 + nY;       
  Y3  = Y2 + nY;       
  Y4  = Y3 + nY;       
  Y5  = Y4 + nY;       
  X5f = Xp + nX;       
  X3f = X5f - 2;         
  X1f = X3f - 2;         
  
  for ( ; Xp < X1f; Xp++) {
    for (Xq = Xp + 1; Xq <= X3f; Xq++) {
      for (Xr = Xq + 1; Xr < X3f; Xr++) {
        for (Xs = Xr + 1; Xs <= X5f; Xs++) {
          for (Xt = Xs + 1; Xt < X5f; Xt++) {
            *Y1++ = *Xp;  
            *Y2++ = *Xq;
            *Y3++ = *Xr;
            *Y4++ = *Xs;
            *Y5++ = *Xt;  
          }
        }
      }
    }
  }
  
  return;
}

void ElemK_8Byte(double *Xp, mwSize nX, double *Yp, mwSize nY, mwSize K)
{
  double *Xq, *Xr,
         **A, **Ap, **Af, *Yq,
         **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
    
  if ((B = (double **) mxMalloc(K * sizeof(double *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (double **) mxMalloc(K * sizeof(double *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  if ((A = (double **) mxMalloc(K * sizeof(double *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;         
       }
       (*Bq)++;           
    }
    
    do {
       Bq--;
       Fp--;
       if (Bq < B) {   
          mxFree(A); 
          mxFree(B);
          mxFree(F);
          return;        
       }
    } while (*Bq >= *Fp); 
    
    (*Bq)++;                
    while (++Bq < Bf) {   
      *Bq = *(Bq - 1) + 1; 
    }
  } 
  
}

void ElemK_4Byte(INT32_T *Xp, mwSize nX, INT32_T *Yp, mwSize nY, mwSize K)
{
  INT32_T *Xq, *Xr,
          **A, **Ap, **Af, *Yq,
          **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
    
  if ((B = (INT32_T **) mxMalloc(K * sizeof(INT32_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (INT32_T **) mxMalloc(K * sizeof(INT32_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  if ((A = (INT32_T **) mxMalloc(K * sizeof(INT32_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;           
       }
       (*Bq)++;             
    }
    
    
    do {
       Bq--;
       Fp--;
       if (Bq < B) {        
          mxFree(A);      
          mxFree(B);
          mxFree(F);
          return;           
       }
    } while (*Bq >= *Fp);  
    
    (*Bq)++;                
    while (++Bq < Bf) { 
      *Bq = *(Bq - 1) + 1; 
    }
  }  
  
}

void ElemK_2Byte(INT16_T *Xp, mwSize nX, INT16_T *Yp, mwSize nY, mwSize K)
{
  INT16_T *Xq, *Xr,
          **A, **Ap, **Af, *Yq,
          **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
  
  if ((B = (INT16_T **) mxMalloc(K * sizeof(INT16_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (INT16_T **) mxMalloc(K * sizeof(INT16_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  if ((A = (INT16_T **) mxMalloc(K * sizeof(INT16_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;          
       }
       (*Bq)++;            
    }
    
    do {
       Bq--;
       Fp--;
       if (Bq < B) {       
          mxFree(A);     
          mxFree(B);
          mxFree(F);
          return;           
       }
    } while (*Bq >= *Fp);  
    
    (*Bq)++;                
    while (++Bq < Bf) {   
      *Bq = *(Bq - 1) + 1;
    }
  }
  
}

void ElemK_1Byte(INT8_T *Xp, mwSize nX, INT8_T *Yp, mwSize nY, mwSize K)
{
  INT8_T  *Xq, *Xr,
          **A, **Ap, **Af, *Yq,
          **B, **Bp, **Bf, **Bq, **lastB, **F, **Fp, **lastF;
    
  if ((B = (INT8_T **) mxMalloc(K * sizeof(INT8_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [B].");
  }
  if ((F = (INT8_T **) mxMalloc(K * sizeof(INT8_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [F].");
  }

  Bf = B + K;
  Fp = F;
  Xq = Xp;
  Xr = Xp + nX - K + 1;
  for (Bp = B; Bp < Bf; ) {
     *Bp++ = Xq++;
     *Fp++ = Xr++;
  }
  
  if ((A = (INT8_T **) mxMalloc(K * sizeof(INT8_T *))) == NULL) {
     mexErrMsgIdAndTxt("JSimon:VChooseK:NoMemory",
                       "Cannot get memory for [A].");
  }
  
  Af = A + K;
  Yq = Yp;
  for (Ap = A; Ap < Af; Ap++) {
     *Ap = Yq;
     Yq += nY;
  }
  
  lastB = B + K - 1;
  lastF = F + K - 1;
  while (1) {
    Bq = lastB;
    Fp = lastF;
    while (*Bq < *Fp) {
       for (Bp = B, Ap = A; Ap < Af; Ap++) {
         **Ap = **Bp++;
         (*Ap)++;          
       }
       (*Bq)++;             
    }
    
    do {
       Bq--;
       Fp--;
       if (Bq < B) {        
          mxFree(A);      
          mxFree(B);
          mxFree(F);
          return;           
       }
    } while (*Bq >= *Fp); 
    
    (*Bq)++;                
    while (++Bq < Bf) { 
      *Bq = *(Bq - 1) + 1; 
    }
  }  
  
}
