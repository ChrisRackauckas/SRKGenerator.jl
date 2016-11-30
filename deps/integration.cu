// filename: integration.cu
// Performs the inner integration loop
extern "C"
{
  __device__ int innerFunc(const float *coefs,const float i,const float isq2,const float isq3,const float isq4,const float isq5,const float isq6,const float isq7,const float isq8,const float j,const float jsq2,const float jsq3,const float jsq4,const float jsq5,const float jsq6,const float jsq7,const float jsq8)
  {
    return abs(coefs[0]*(jsq2) + coefs[1]*(jsq3) + coefs[2]*(jsq4) + coefs[3]*(jsq5) + coefs[4]*jsq6 + coefs[5]*jsq7 + coefs[6]*jsq8 + coefs[7]*(i) + coefs[8]*(i)*(jsq2) + coefs[9]*i*jsq3 + coefs[10]*(i)*(jsq4) + coefs[11]*i*jsq5 + coefs[12]*(i)*(jsq6) + coefs[13]*i*jsq7 + coefs[14]*(isq2) + coefs[15]*(isq2)*(jsq2) + coefs[16]*isq2*jsq3 + coefs[17]*(isq2)*(jsq4) + coefs[18]*isq2*jsq5 + coefs[19]*(isq2)*(jsq6) + coefs[20]*(isq3) + coefs[21]*(isq3)*(jsq2) + coefs[22]*isq3*jsq3 + coefs[23]*(isq3)*(jsq4) + coefs[24]*isq3*jsq5 + coefs[25]*(isq4) + coefs[26]*(isq4)*(jsq2) + coefs[27]*isq4*jsq3 + coefs[28]*(isq4)*(jsq4) + coefs[29]*(isq5) + coefs[30]*(isq5)*(jsq2) + coefs[31]*isq5*jsq3+ coefs[32]*(isq6) + coefs[33]*(isq6)*(jsq2) + coefs[34]*(isq7) + coefs[35]*(isq8))<1;
  }
    __global__ void integration(const float *coefs, const float *iArr, const float *jArr, const int sizei, const int sizej, const int equalDiv,const int startIdx, int *tmp)
    {
        int index = threadIdx.x + blockIdx.x * blockDim.x;
        int globalIndex = index*equalDiv+startIdx;
        int loopInd;
        float i;
        float j;
        float isq2;
        float isq3;
        float isq4;
        float isq5;
        float isq6;
        float isq7;
        float isq8;
        float jsq2;
        float jsq3;
        float jsq4;
        float jsq5;
        float jsq6;
        float jsq7;
        float jsq8;
        int ans = 0;
        for(loopInd=0;loopInd<equalDiv;loopInd=loopInd+1){
          i = iArr[(globalIndex+loopInd)/sizej];
          j = jArr[(globalIndex+loopInd)%sizej];
          if(globalIndex+loopInd >= sizei*sizej){
            break;
          }
          if((globalIndex+loopInd)%sizej==0 || loopInd==0){
            isq2 = i*i;
            isq3 = i*isq2;
            isq4 = isq2*isq2;
            isq5 = i*isq4;
            isq6 = isq4*isq2;
            isq7 = i*isq6;
            isq8 = isq4*isq4;
          }
          jsq2 = j*j;
          jsq3 = j*jsq2;
          jsq4 = jsq2*jsq2;
          jsq5 = j*jsq4;
          jsq6 = jsq2*jsq4;
          jsq7 = j*jsq6;
          jsq8 = jsq4*jsq4;
          ans = ans + innerFunc(coefs,i,isq2,isq3,isq4,isq5,isq6,isq7,isq8,j,jsq2,jsq3,jsq4,jsq5,jsq6,jsq7,jsq8);
        }
        tmp[index] = ans;
    }
}
