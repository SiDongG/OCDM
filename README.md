# OCDM
MATLAB implementation of Orthogonal Chirp Division Multiplexing w Discrete Fresnel Transform 


### OCDM folder: 

OCDM.mat generates 4,16,64-QAM OCDM transmission bit rate under ZF and MMSE equalizations

Transmitter.m, Channel.m, Receiver.m are Transmitter, multipath channel, receiver architecture simulations respectively 

Implementation algorithm based on Seminal Paper on OCDM:https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7523229

--Remark1: Equation 34, Instead of Complex Conjugate of diagonal CFR matrix, Conjugate Transpose should be used fpr consistency, although result does not change  

--Remark2: Currently Only support even number of sub-carriers, odd number Fresnel transform not implemented (Update: fixed, but equation 30 odd number Zadoff-Chu Sequence is changed to k(k+1) for it to work, needs mathematically proof of relation under equation 29) 

--Remark3: Equation 30 is questionable, instead of k, use (N-k+1) as index (i.e.Reverse the diagonal)


### DFT-P-OFDM folder: 

DFTPOFDM.mat generates single carrier DFT precoded transmission and calculates PAPR before and after DFT-Precoding. around 1.5dB decrease. 

Algorithm based on WINNER project DFT-P-OFDM:https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5357586&tag=1

--Remark1: Did not implement entire transmission chain and plot BER, the principle is simple SCM transmission


### OCDM-CFO folder:

Simulate performance of OCDM under subcarrier frequency offset and blind CSI w ZF and MMSE equalization 

CFO mitigation and channel estimation implementation: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9638964

CFOOCDM.m demonstrates successful CFO compensation and channel estimation with OCDM transmission, Transceiver architecture in OCDMTxRx.m

--Remark1: The model assumes time-invariant channel, but principle can be applied to quasi-static channel easily

--Remark2: Equation 4, Progressive CFO based on Block number term is omitted, block index dropped during equalization

--Remark3: Equation 14, definition of the matrix has an error on H12, should be (K-Lh-1) zeros

--Remark4: Pilots are longer than necessary, but since zero-paddling pilots are used, it is not energy-inefficient 

### MIMO-OCDM
 2 by 2 MIMO-OCDM linear equalization 
 
### FDPrecodedMIMO-OCDM
 2 by 2 Frequency domain precoded MIMO-OCDM, ideas modified based on and extended from: https://ieeexplore.ieee.org/document/1296644
