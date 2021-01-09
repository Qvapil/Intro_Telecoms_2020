import matplotlib.pyplot as plt
import numpy as np
import math
import graycode
from scipy.io import wavfile

#constants
A=1
AMi=2
fmi=5

#1.a(i): sample signal y(t) with sampling frequency fs1=20fm
def yi_func(t):
    return A*np.cos(2*math.pi*fmi*t)*np.cos(2*math.pi*(AMi+2)*fmi*t)

fs1i=20*fmi
t1i=np.linspace(0,4/fmi,int(4*fs1i/fmi))
y1i=yi_func(t1i)

fig1i, ax1i = plt.subplots()
plt.stem(t1i,y1i,markerfmt='.',basefmt='-k')
ax1i.set_title('y(t) sampled with frequency $fs_1$=20$f_m$')
ax1i.set_xlabel('Time (ms)')
ax1i.set_ylabel('Amplitude (V)')

#1.a(ii): sample signal y(t) with sampling frequency fs2=100fm
fs2i=100*fmi
t2i=np.linspace(0,4/fmi,int(4*fs2i/fmi))
y2i=yi_func(t2i)

fig2i, ax2i = plt.subplots()
plt.stem(t2i,y2i,markerfmt='.',basefmt='-k')
ax2i.set_title('y(t) sampled with frequency $fs_2$=100$f_m$')
ax2i.set_xlabel('Time (ms)')
ax2i.set_ylabel('Amplitude (V)')

#1.a(iii): plot both samples
fig3i, ax3i=plt.subplots()
plt.stem(t2i,y2i,linefmt='-r',markerfmt='.r',basefmt='-k',label='$fs_2$=100$f_m$')
plt.stem(t1i,y1i,linefmt='-b',markerfmt='.b',basefmt='-k',label='$fs_1$=20$f_m$')
ax3i.set_title('y(t) sampled with frequencies $fs_1$ and $fs_2$')
ax3i.set_xlabel('Time (ms)')
ax3i.set_ylabel('Amplitude (V)')
ax3i.legend()

#1.b: sample y(t) with sampling frequency fs3=5fm
fs3i=5*fmi
t3i=np.linspace(0,4/fmi,int(4*fs3i/fmi))
y3i=yi_func(t3i)

fig4i, ax4i = plt.subplots()
plt.stem(t3i,y3i,markerfmt='.',basefmt='-k')
ax4i.set_title('y(t) sampled with frequency $fs_3$=5$f_m$')
ax4i.set_xlabel('Time (ms)')
ax4i.set_ylabel('Amplitude (V)')

#2.a: quantizing y(t) sampled with sampling frequency fs1=20fm

Ri=5
Li=2**Ri
Di=2*(max(y1i))/(Li-1)
y1i_quant=Di*(np.floor(y1i/Di)+1/2)

gray_code=['00000','00001','00011','00010','00110','00111','00101','00100','01100','01101','01111','01110','01010','01011','01001','01000','11000','11001','11011','11010','11110','11111','11101','11100','10100','10101','10111','10110','10010','10011','10001','10000']
quant2levelsi=np.linspace(-1,1,32)
fig5i, ax5i = plt.subplots()
plt.stem(t1i,y1i_quant,markerfmt='.',basefmt='-k')
plt.grid(axis='y')
plt.yticks(quant2levelsi,gray_code)
ax5i.set_title('Quantized y(t) sampled with frequency $fs_1$=20$f_m$')
ax5i.set_xlabel('Time (ms)')
ax5i.set_ylabel('Gray Code')

#2.b(i): standard deviation of error of first 10 samples
errori=y1i-y1i_quant
print('Standard Deviation of error of first 10 samples: {}'.format(np.std(errori[:10])))

#2.b(ii): standard deviation of error of first 20 samples
print('Standard Deviation of error of first 20 samples: {}'.format(np.std(errori[:20])))

#2.b(iii): SNR of (i), (ii) and theoretical value
Fi10=max(errori[:10])/np.sqrt(np.mean(errori[:10]**2))
print('SNR of first 10 samples: {}'.format(6.02*Ri-20*np.log10(Fi10)+4.77))

Fi20=max(errori[:20])/np.sqrt(np.mean(errori[:20]**2))
print('SNR of first 20 samples: {}'.format(6.02*Ri-20*np.log10(Fi20)+4.77))

print('Theoretical SNR: {}'.format(20*np.log10(2**Ri)))

#2.c: polar nrz of 1 period of quantized y(t)
y1i_quant_1period=y1i_quant[:int(len(y1i_quant)/4)]

# initially we tried this:  dicti=dict(zip(np.linspace(-1,1,32),gray_code))
# but python doesnt deal with the precise floats well
# and didnt recognise some values of y as equal to the levels of the linspace
# so we had to manually set the gray levels to the y values
gray_code_2i=['00000','00111','00101','00100','01101','01011','01001','01000','11000','11001','11110','11111','10100','10101','10111','10000']
dicti=dict(zip(sorted(set(y1i_quant_1period)),gray_code_2i))

bitsnrzi=''
for i in range(len(y1i_quant_1period)):
     bitsnrzi+=dicti.get(y1i_quant_1period[i])

def modulatehli(bits,a):
    modsig=np.zeros(len(bits))
    for i in range(len(bits)):
        if bits[i]==0 or bits[i]=='0':
            modsig[i]=-a
        elif bits[i]==1 or bits[i]=='1':
            modsig[i]=a
    return modsig

polar_nrz_i=modulatehli(bitsnrzi,fmi)

fig6i, ax6i = plt.subplots()
plt.plot(polar_nrz_i,drawstyle='steps-post')
ax6i.set_title('Polar NRZ (1 period of quantized y(t))')
ax6i.set_xlabel('Time (ms)')
ax6i.set_ylabel('Amplitude(V)')

#3a:B-PAM of random bitstream
bitstream2i=np.random.randint(2,size=46)
print('Bits for ex 3,4:')
print(bitstream2i)
A3i=5
Tb=0.5
t3ai=np.linspace(0,46*Tb,46)
Ei=A3i**2*Tb

BPAMi=modulatehli(bitstream2i,A3i)

def energyi(a):
    return np.sign(a)*np.sqrt(a**2*Tb)

fig7i, ax7i = plt.subplots()
plt.plot(t3ai,BPAMi,drawstyle='steps-post')
ax7i.set_title('BPAM of randomised bitstream')
ax7i.set_xlabel('Time (s)')
ax7i.set_ylabel('Amplitude(V)')

#3b: constellation of no-noise bitstream
fig8i, ax8i = plt.subplots()
plt.plot(energyi(BPAMi),np.imag(BPAMi),'*r')
plt.xlim(-7,7)
plt.ylim(-7,7)
ax8i.set_aspect('equal')
ax8i.set_title('Constellation Diagram of BPAM')
ax8i.set_xlabel('In phase')
ax8i.set_ylabel('Quadrature')

#3c: adding noise to BPAM signal
def awgn_ampli(length,SNR,E):
    N0=E/20**(SNR/10)
    noise=np.random.normal(0,np.sqrt(N0),length)
    return noise

noisyBPAMi1=BPAMi+awgn_ampli(len(BPAMi),5,Ei)

fig9i, ax9i = plt.subplots()
plt.plot(t3ai,noisyBPAMi1,drawstyle='steps-post')
ax9i.set_title('Noisy BPAM signal with SNR=5db')
ax9i.set_xlabel('Time (s)')
ax9i.set_ylabel('Amplitude (V)')

noisyBPAMi2=BPAMi+awgn_ampli(len(BPAMi),15,Ei)

fig10i, ax10i = plt.subplots()
plt.plot(t3ai,noisyBPAMi2,drawstyle='steps-post')
ax10i.set_title('Noisy BPAM signal with SNR=15db')
ax10i.set_xlabel('Time (s)')
ax10i.set_ylabel('Amplitude (V)')

#3d: constellation of noisy signals
fig11i, ax11i = plt.subplots()
plt.plot(energyi(noisyBPAMi1),energyi(awgn_ampli(len(BPAMi),5,Ei)),'o',label='signal with noise')
plt.plot(energyi(BPAMi),np.imag(BPAMi),'*r',label='signal without noise')
plt.xlim(-7,7)
plt.ylim(-7,7)
ax11i.set_aspect('equal')
ax11i.set_title('Constellation Diagram of noisy BPAM signal with SNR=5db')
ax11i.set_xlabel('In phase')
ax11i.set_ylabel('Quadrature')
ax11i.legend()

fig12i, ax12i = plt.subplots()
plt.plot(energyi(noisyBPAMi2),energyi(awgn_ampli(len(BPAMi),15,Ei)),'o',label='signal with noise')
plt.plot(energyi(BPAMi),np.imag(BPAMi),'*r',label='signal without noise')
plt.xlim(-7,7)
plt.ylim(-7,7)
ax12i.set_aspect('equal')
ax12i.set_title('Constellation Diagram of noisy BPAM signal with SNR=15db')
ax12i.set_xlabel('In phase')
ax12i.set_ylabel('Quadrature')
ax12i.legend()

#3e: bit error rate of BPAM
lengthi=100000
bitstream3i=np.random.randint(2,size=lengthi)

BPAMei=modulatehli(bitstream3i,A3i)

SNRei=np.arange(0,16)
N0ei=Ei/20**(SNRei/10)
noisysignal3i=np.zeros((16,lengthi))
errors3i=np.zeros(16)
BER3i=np.zeros(16)

def bpamErrorsi(signal,noisysignal):
    errors=0
    for i in range(len(signal)):
        if noisysignal[i]*signal[i]<0:
            errors+=1
    return errors

for i in range(16):
    noisysignal3i[i]=BPAMei+awgn_ampli(len(BPAMei),i,Ei)
    errors3i[i]=bpamErrorsi(BPAMei,noisysignal3i[i])
    if errors3i[i]!=0:
        BER3i[i]=np.log10(errors3i[i]/lengthi)
    else:
        BER3i[i]=-50

Q3i=np.zeros(16)
for i in SNRei:
    for j in np.linspace(np.sqrt(2*Ei/N0ei[i]),1000,1000):
        Q3i[i]+=math.exp(-j**2/2)/np.sqrt(2*math.pi)

fig13i, ax13i = plt.subplots()
plt.plot(BER3i,'or',label='Empirical')
plt.plot(np.log10(Q3i),label='Theoretical')
plt.ylim(bottom=-40)
ax13i.set_title('BER of BPAM')
ax13i.set_xlabel('SNR (db)')
ax13i.set_ylabel('$log_{10}$BER')
ax13i.legend()

#4a: constellation of QPSK
QPSKxi=np.zeros(23)
QPSKyi=np.zeros(23)
a4i=np.sqrt(Ei)/np.sqrt(2)
for i in range(0,46,2):
    if bitstream2i[i]==0 and bitstream2i[i+1]==0:
        QPSKxi[int(i/2)]=a4i
        QPSKyi[int(i/2)]=a4i
    elif bitstream2i[i]==0 and bitstream2i[i+1]==1:
        QPSKxi[int(i/2)]=-a4i
        QPSKyi[int(i/2)]=a4i
    elif bitstream2i[i]==1 and bitstream2i[i+1]==1:
        QPSKxi[int(i/2)]=-a4i
        QPSKyi[int(i/2)]=-a4i
    elif bitstream2i[i]==1 and bitstream2i[i+1]==0:
        QPSKxi[int(i/2)]=a4i
        QPSKyi[int(i/2)]=-a4i


fig13i, ax13i = plt.subplots()
plt.plot(QPSKxi,QPSKyi,'*r')
plt.text(2.6,2.6,'00',fontsize=9)
plt.text(-3,2.6,'01',fontsize=9)
plt.text(-3,-3,'11',fontsize=9)
plt.text(2.6,-3,'10',fontsize=9)
plt.xlim(-6,6)
plt.ylim(-6,6)
ax13i.set_aspect('equal')
ax13i.set_title('Constellation Diagram of QPSK modulated signal')
ax13i.set_xlabel('In phase')
ax13i.set_ylabel('Quadrature')

#4b: adding noise to QPSK signal and plot constellation diagrams
noisyQPSKxi1=QPSKxi+energyi(awgn_ampli(len(QPSKxi),5,Ei))
noisyQPSKyi1=QPSKyi+energyi(awgn_ampli(len(QPSKyi),5,Ei))

noisyQPSKxi2=QPSKxi+energyi(awgn_ampli(len(QPSKxi),15,Ei))
noisyQPSKyi2=QPSKyi+energyi(awgn_ampli(len(QPSKyi),15,Ei))

fig14i, ax14i = plt.subplots()
plt.scatter(noisyQPSKxi1,noisyQPSKyi1,label='signal with noise')
plt.plot(QPSKxi,QPSKyi,'*r',label='signal without noise')
plt.xlim(-6,6)
plt.ylim(-6,6)
ax14i.set_aspect('equal')
ax14i.set_title('Constellation Diagram of noisy QPSK signal with SNR=5db')
ax14i.set_xlabel('In phase')
ax14i.set_ylabel('Quadrature')
ax14i.legend()

fig15i, ax15i = plt.subplots()
plt.scatter(noisyQPSKxi2,noisyQPSKyi2,label='signal with noise')
plt.plot(QPSKxi,QPSKyi,'*r',label='signal without noise')
plt.xlim(-6,6)
plt.ylim(-6,6)
ax15i.set_aspect('equal')
ax15i.set_title('Constellation Diagram of noisy QPSK signal with SNR=15db')
ax15i.set_xlabel('In phase')
ax15i.set_ylabel('Quadrature')
ax15i.legend()

#4c: bit error rate of QPSK
lengthi=10000
bitstream4i=np.random.randint(2,size=lengthi)
QPSKcxi=np.zeros(int(lengthi/2))
QPSKcyi=np.zeros(int(lengthi/2))

for i in range(0,int(lengthi/2),2):
    if bitstream4i[i]==0 and bitstream4i[i+1]==0:
        QPSKcxi[int(i/2)]=a4i
        QPSKcyi[int(i/2)]=a4i
    elif bitstream4i[i]==0 and bitstream4i[i+1]==1:
        QPSKcxi[int(i/2)]=-a4i
        QPSKcyi[int(i/2)]=a4i
    elif bitstream4i[i]==1 and bitstream4i[i+1]==1:
        QPSKcxi[int(i/2)]=-a4i
        QPSKcyi[int(i/2)]=-a4i
    elif bitstream4i[i]==1 and bitstream4i[i+1]==0:
        QPSKcxi[int(i/2)]=a4i
        QPSKcyi[int(i/2)]=-a4i

SNR4i=np.arange(0,16)
N04i=Ei/20**(SNR4i/10)
noisysignalxi=np.zeros((16,int(lengthi/2)))
noisysignalyi=np.zeros((16,int(lengthi/2)))
errors4i=np.zeros(16)
BER4i=np.zeros(16)

def qpskErrorsi(sigx,sigy,noisyx,noisyy):
    errors=0
    for i in range(len(sigx)):
        if sigx[i]*noisyx[i]<0 or sigy[i]*noisyy[i]<0:
            errors+=1
    return errors

for i in SNR4i:
    noisysignalxi[i]=QPSKcxi+energyi(awgn_ampli(len(QPSKcxi),i,Ei))
    noisysignalyi[i]=QPSKcyi+energyi(awgn_ampli(len(QPSKcyi),i,Ei))
    errors4i[i]=qpskErrorsi(QPSKcxi,QPSKcyi,noisysignalxi[i],noisysignalyi[i])
    if errors4i[i]!=0:
        BER4i[i]=np.log10(errors4i[i]/(lengthi/2))
    else:
        BER4i[i]=-50


Q4=np.zeros(16)
for i in SNR4i:
    for j in np.linspace(np.sqrt(Ei/N04i[i]),1000,1000):
        Q4[i]+=math.exp(-j**2/2)/np.sqrt(2*math.pi)

fig16i, ax16i = plt.subplots()
plt.plot(SNR4i,BER4i,'or',label='Empirical')
plt.plot(SNR4i,np.log10(Q4),label='Theoretical')
plt.ylim(bottom=-20)
ax16i.set_title('BER of QPSK')
ax16i.set_xlabel('SNR (db)')
ax16i.set_ylabel('$log_{10}$BER')
ax16i.legend()


#4.d(i): read file as binary
with open('shannon_odd.txt', 'r') as file:
    datai = file.read().replace('\n', '')

databini=bin(int.from_bytes(datai.encode(),'big'))
databini=databini.replace('b','')
datadeci=np.zeros(int(len(databini)/8))
for i in range(0,len(databini),8):
    datadeci[int(i/8)]=int(databini[i:i+8],2)


#4.d(ii): quantize with 8 bits
def quantize(x,min,max,levels):
    x_normalize=(x-min)*(levels-1)/(max-min)
    x_normalize[x_normalize>levels-1]=levels-1
    x_normalize[x_normalize<0]=0
    x_normalize_quant=np.around(x_normalize)
    x_quant=(x_normalize_quant)*(max-min)/(levels-1)+min
    return x_quant

dataiquant=quantize(datadeci,max(datadeci),min(datadeci),256)
quantlevelsi=np.linspace(min(datadeci),max(datadeci),256)

graycode8=graycode.gen_gray_codes(8)

dict4i=dict(zip(np.linspace(min(datadeci),max(datadeci),256),graycode8))

grayvaluei=np.zeros(len(dataiquant))
for i in range(len(dataiquant)):
    for j in range(256):
        if math.isclose(dataiquant[i],quantlevelsi[j],abs_tol=0.0001):
            grayvaluei[i]=dict4i.get(quantlevelsi[j])

fig17i, ax17i = plt.subplots()
plt.plot(dataiquant,'.')
plt.yticks(np.linspace(min(datadeci),max(datadeci),256),graycode8)
ax17i.grid(axis='y')
ax17i.set_title('Quantized text')
ax17i.set_ylabel('Gray code (in decimal)')

#4d(iii): QPSK modulation
bits4di=''
for i in range(len(grayvaluei)):
    bnr = bin(int(grayvaluei[i])).replace('0b','')
    x = bnr[::-1]
    while len(x) < 8:
        x += '0'
    bnr = x[::-1]
    bits4di+=str(bnr)

textQPSKxi=np.zeros(int(len(bits4di)/2))
textQPSKyi=np.zeros(int(len(bits4di)/2))
a4di=1/np.sqrt(2)

for i in range(0,len(bits4di),2):
    if bits4di[i]=='0' and bits4di[i+1]=='0':
        textQPSKxi[int(i/2)]=a4di
        textQPSKyi[int(i/2)]=a4di
    elif bits4di[i]=='0' and bits4di[i+1]=='1':
        textQPSKxi[int(i/2)]=-a4di
        textQPSKyi[int(i/2)]=a4di
    elif bits4di[i]=='1' and bits4di[i+1]=='1':
        textQPSKxi[int(i/2)]=-a4di
        textQPSKyi[int(i/2)]=-a4di
    elif bits4di[i]=='1' and bits4di[i+1]=='0':
        textQPSKxi[int(i/2)]=a4di
        textQPSKyi[int(i/2)]=-a4di

fig18i, ax18i = plt.subplots()
plt.plot(textQPSKxi,textQPSKyi,'*r')
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
ax18i.set_aspect('equal')
ax18i.set_title('Constellation Diagram of QPSK modulated text')
ax18i.set_xlabel('In phase')
ax18i.set_ylabel('Quadrature')

#4.d.(iv):adding noise
Tb=1
Eti=1

noisytextQPSKxi1=textQPSKxi+energyi(awgn_ampli(len(textQPSKxi),5,Eti))
noisytextQPSKyi1=textQPSKyi++energyi(awgn_ampli(len(textQPSKyi),5,Eti))

noisytextQPSKxi2=textQPSKxi+energyi(awgn_ampli(len(textQPSKxi),15,Eti))
noisytextQPSKyi2=textQPSKyi+energyi(awgn_ampli(len(textQPSKyi),15,Eti))

#4.d.(v):constellations of noisy signals
fig19i, ax19i = plt.subplots()
plt.scatter(noisytextQPSKxi1,noisytextQPSKyi1,label='signal with noise')
plt.plot(textQPSKxi,textQPSKyi,'*r',label='signal without noise')
plt.xlim(-2,2)
plt.ylim(-2,2)
ax19i.set_aspect('equal')
ax19i.set_title('Constellation Diagram of noisy QPSK text with SNR=5db')
ax19i.set_xlabel('In phase')
ax19i.set_ylabel('Quadrature')
ax19i.legend()

fig20i, ax20i = plt.subplots()
plt.scatter(noisytextQPSKxi2,noisytextQPSKyi2,label='signal with noise')
plt.plot(textQPSKxi,textQPSKyi,'*r',label='signal without noise')
plt.xlim(-2,2)
plt.ylim(-2,2)
ax20i.set_aspect('equal')
ax20i.set_title('Constellation Diagram of noisy QPSK text with SNR=15db')
ax20i.set_xlabel('In phase')
ax20i.set_ylabel('Quadrature')
ax20i.legend()

# 4.d.(vi): BER
errors4di=np.zeros(2)
for i in range(len(textQPSKxi)):
    if noisytextQPSKxi1[i]*textQPSKxi[i]<0 or noisytextQPSKyi1[i]*textQPSKyi[i]<0:
        errors4di[0]+=1
    if noisytextQPSKxi2[i]*textQPSKxi[i]<0 or noisytextQPSKyi2[i]*textQPSKyi[i]<0:
        errors4di[1]+=1

Qt=np.zeros(2)
N0ti1=Eti/20**(5/10)
N0ti2=Eti/20**(15/10)
for i in np.linspace(np.sqrt(Eti/N0ti1),1000,1000):
    Qt[0]+=math.exp(-i**2/2)/np.sqrt(2*math.pi)
for i in np.linspace(np.sqrt(Eti/N0ti2),1000,1000):
    Qt[1]+=math.exp(-i**2/2)/np.sqrt(2*math.pi)

print('Experimental BER for 5db noise: {}'.format(errors4di[0]/len(textQPSKxi)))
print('Experimental BER for 15db noise: {}'.format(errors4di[1]/len(textQPSKxi)))
print('Theoretical BER for 5db noise: {}'.format(Qt[0]))
print('Theoretical BER for 15db noise: {}'.format(Qt[1]))

#4.d.(vii): remake text
def qpskDemodulate(x,y):
    bits=''
    for i in range(len(x)):
        if x[i]>0 and y[i]>0:
            bits+='00'
        if x[i]<0 and y[i]>0:
            bits+='01'
        if x[i]<0 and y[i]<0:
            bits+='11'
        if x[i]>0 and y[i]<0:
            bits+='10'
    return bits

newbits5i=qpskDemodulate(noisytextQPSKxi1,noisytextQPSKyi1)
newbits15i=qpskDemodulate(noisytextQPSKxi2,noisytextQPSKyi2)

newgrayvalue5i=np.zeros(int(len(newbits5i)/8))
newgrayvalue15i=np.zeros(int(len(newbits5i)/8))

for i in range(0,len(newbits5i),8):
    newgrayvalue5i[int(i/8)]=int(newbits5i[i:i+8],2)
    newgrayvalue15i[int(i/8)]=int(newbits15i[i:i+8],2)

newdec5i=np.zeros(len(newgrayvalue5i))
newdec15i=np.zeros(len(newgrayvalue15i))
inv_dict4i = {v: k for k, v in dict4i.items()}
for i in range(len(newgrayvalue5i)):
    newdec5i[i]=np.round(inv_dict4i.get(newgrayvalue5i[i]))
    newdec15i[i]=np.round(inv_dict4i.get(newgrayvalue15i[i]))

text5i=''
text15i=''
for i in range(len(newgrayvalue5i)):
    text5i+=chr(int(newdec5i[i]))
    text15i+=chr(int(newdec15i[i]))

f = open("shannon_odd_5dbi.txt", "w")
f.write(text5i)
f.close()

f = open("shannon_odd_15dbi.txt", "w")
f.write(text15i)
f.close()

#5a: plot sound file
sampleratei,data5i=wavfile.read('soundfile1_lab2.wav')

fig21i, ax21i = plt.subplots()
plt.plot(data5i)
ax21i.grid(axis='y')
ax21i.set_title('Sound signal')
ax21i.set_xlabel('Samples')

#5b: quantize sound signal
Ri=8
Li=2**Ri
maxi=max(data5i)
Di=2*maxi/(Li-1)
data5iquant=Di*(np.floor(data5i/Di)+1/2)

quantlevels5i=np.linspace(-max(data5i),max(data5i),256)
graycode8=graycode.gen_gray_codes(8)

fig22i, ax22i = plt.subplots()
plt.plot(data5iquant,'.')
plt.yticks(quantlevels5i,graycode8)
ax22i.grid(axis='y')
ax22i.set_title('Quantized sound signal')
ax22i.set_ylabel('Gray code (in decimal)')

#5c: qpsk modulation of sound signal
grayvalues5i=np.zeros(len(data5iquant))
for i in range(len(data5iquant)):
    grayvalues5i[i]=graycode8[int((data5iquant[i]+maxi)/Di)]

bits5i=''
for i in range(len(grayvalues5i)):
    bnr = bin(int(grayvalues5i[i])).replace('0b','')
    x = bnr[::-1]
    while len(x) < 8:
        x += '0'
    bnr = x[::-1]
    bits5i+=str(bnr)

soundQPSKxi=np.zeros(int(len(bits5i)/2))
soundQPSKyi=np.zeros(int(len(bits5i)/2))
a5i=1/np.sqrt(2)

for i in range(0,len(bits5i),2):
    if bits5i[i]=='0' and bits5i[i+1]=='0':
        soundQPSKxi[int(i/2)]=a5i
        soundQPSKyi[int(i/2)]=a5i
    elif bits5i[i]=='0' and bits5i[i+1]=='1':
        soundQPSKxi[int(i/2)]=-a5i
        soundQPSKyi[int(i/2)]=a5i
    elif bits5i[i]=='1' and bits5i[i+1]=='1':
        soundQPSKxi[int(i/2)]=-a5i
        soundQPSKyi[int(i/2)]=-a5i
    elif bits5i[i]=='1' and bits5i[i+1]=='0':
        soundQPSKxi[int(i/2)]=a5i
        soundQPSKyi[int(i/2)]=-a5i

fig23i, ax23i = plt.subplots()
plt.plot(soundQPSKxi,soundQPSKyi,'*r')
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
ax23i.set_aspect('equal')
ax23i.set_title('Constellation Diagram of QPSK modulated sound')
ax23i.set_xlabel('In phase')
ax23i.set_ylabel('Quadrature')

#5d: adding noise
symbolxi=[a5i,-a5i,a5i,-a5i]
symbolyi=[a5i,a5i,-a5i,-a5i]

Tb=1
Eti=1

noisysoundQPSKxi1=soundQPSKxi+energyi(awgn_ampli(len(soundQPSKxi),4,Eti))
noisysoundQPSKyi1=soundQPSKyi++energyi(awgn_ampli(len(soundQPSKyi),4,Eti))

noisysoundQPSKxi2=soundQPSKxi+energyi(awgn_ampli(len(soundQPSKxi),14,Eti))
noisysoundQPSKyi2=soundQPSKyi+energyi(awgn_ampli(len(soundQPSKyi),14,Eti))

#5e: constellations of noisy signals
fig24i, ax24i = plt.subplots()
plt.scatter(noisysoundQPSKxi1,noisysoundQPSKyi1,marker='.',label='signal with noise')
plt.plot(symbolxi,symbolyi,'*r',label='signal without noise')
plt.xlim(-3,3)
plt.ylim(-3,3)
ax24i.set_aspect('equal')
ax24i.set_title('Constellation Diagram of noisy QPSK sound with SNR=4db')
ax24i.set_xlabel('In phase')
ax24i.set_ylabel('Quadrature')
ax24i.legend()

fig25i, ax25i = plt.subplots()
plt.scatter(noisysoundQPSKxi2,noisysoundQPSKyi2,marker='.',label='signal with noise')
plt.plot(symbolxi,symbolyi,'*r',label='signal without noise')
plt.xlim(-2,2)
plt.ylim(-2,2)
ax25i.set_aspect('equal')
ax25i.set_title('Constellation Diagram of noisy QPSK sound with SNR=14db')
ax25i.set_xlabel('In phase')
ax25i.set_ylabel('Quadrature')
ax25i.legend()

#5st: BER of noisy sound
errors5i=np.zeros(2)
for i in range(len(soundQPSKxi)):
    if noisysoundQPSKxi1[i]*soundQPSKxi[i]<0 or noisysoundQPSKyi1[i]*soundQPSKyi[i]<0:
        errors5i[0]+=1
    if noisysoundQPSKxi2[i]*soundQPSKxi[i]<0 or noisysoundQPSKyi2[i]*soundQPSKyi[i]<0:
        errors5i[1]+=1

Qs=np.zeros(2)
N0ti1=Eti/20**(4/10)
N0ti2=Eti/20**(14/10)
for i in np.linspace(np.sqrt(Eti/N0ti1),1000,1000):
    Qs[0]+=math.exp(-i**2/2)/np.sqrt(2*math.pi)
for i in np.linspace(np.sqrt(Eti/N0ti2),1000,1000):
    Qs[1]+=math.exp(-i**2/2)/np.sqrt(2*math.pi)

print('Experimental BER for 4db noise: {}'.format(errors5i[0]/len(soundQPSKxi)))
print('Experimental BER for 14db noise: {}'.format(errors5i[1]/len(soundQPSKxi)))
print('Theoretical BER for 4db noise: {}'.format(Qs[0]))
print('Theoretical BER for 14db noise: {}'.format(Qs[1]))

#5z: remaking audio signal
newbits4si=qpskDemodulate(noisysoundQPSKxi1,noisysoundQPSKyi1)
newbits14si=qpskDemodulate(noisysoundQPSKxi2,noisysoundQPSKyi2)

newgrayvalue4si=np.zeros(int(len(newbits4si)/8))
newgrayvalue14si=np.zeros(int(len(newbits4si)/8))

for i in range(0,len(newbits4si),8):
    newgrayvalue4si[int(i/8)]=int(newbits4si[i:i+8],2)
    newgrayvalue14si[int(i/8)]=int(newbits14si[i:i+8],2)

dict5i=dict(zip(quantlevels5i,graycode8))

newdec4si=np.zeros(len(newgrayvalue4si))
newdec14si=np.zeros(len(newgrayvalue14si))
inv_dict5i = {v: k for k, v in dict5i.items()}
for i in range(len(newgrayvalue4si)):
    newdec4si[i]=np.round(inv_dict5i.get(newgrayvalue4si[i]))
    newdec14si[i]=np.round(inv_dict5i.get(newgrayvalue14si[i]))

wavfile.write('soundfile1_lab2_4dbi.wav',sampleratei,newdec4si)
wavfile.write('soundfile1_lab2_14dbi.wav',sampleratei,newdec14si)


plt.show()
