# OCDM
MATLAB implementation of Orthogonal Chirp Division Multiplexing w Discrete Fresnel Transform 


### OCDM folder: 

OCDM.mat generates 4,16,64-QAM OCDM transmission bit rate under ZF and MMSE equalizations

Transmitter.m, Channel.m, Receiver.m are Transmitter, multipath channel, receiver architecture simulations respectively 

Implementation algorithm based on Seminal Paper on OCDM:https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7523229

--Remark1: Equation 34 is erraneous, Instead of Complex Conjugate of diagonal CFR matrix, Conjugate Transpose should be used 

--Remark2: Currently Only support even number of sub-carriers, odd number Fresnel transform not implemented (Update: fixed, but equation 30 odd number Zadoff-Chu Sequence is changed to k(k+1) for it to work, needs mathematically proof of relation under equation 29) 

--Remark3: Equation 30 is questionable, instead of k, use (N-k+1) as index (i.e.Reverse the diagonal)


### DFT-P-OFDM folder: 

DFTPOFDM.mat generates single carrier DFT precoded transmission and calculates PAPR
