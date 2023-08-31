# FFT-msvc
Fast Fourier Transform in C

## Implemented Algorithms
* One-dimensional Discrete Fourier Transform
* One-dimensional Fast Fourier Transform
  * 8-points Decimation-in-time Algorithm
  * Radix-2 Cooley–Tukey Algorithm
  * Radar's Algorithm
  * Hybrid Radix Cooley–Tukey Algorithm

## Algorithms Overview
### One-dimensional Discrete Fourier Transform
> #### Discrete Fourier Transform Algorithm - \$\mathcal{O}(n^2)\$
> $$ X_k = \sum_{n=0}^{N-1}x_n\ e^{-\frac{i2\pi}{N}nk} \qquad k = 0, \cdots, N-1 $$
> 
> **Time Complexity :**\
> n summation to get a frequency sample, a total of n frequency sample \$\rightarrow\$ \$n\times n\$

### One-dimensional Fast Fourier Transform
> #### 8-points Decimate-in-time Algorithm - \$\mathcal{O}(n)\$
>> $$ W_N^0 = e^{-\frac{i2\pi}{N}0} = 1 $$
>>
>> <picture>
>>   <source media="(prefers-color-scheme: dark)" srcset="https://github.com/Yuyu378/FFT-msvc/blob/main/img/8-points%20Decimate-in-time%20Algorithm_black.jpg">
>>   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Yuyu378/FFT-msvc/blob/main/img/8-points%20Decimate-in-time%20Algorithm.jpg">
>>   <img alt="Figure 1" src="https://github.com/Yuyu378/FFT-msvc/blob/main/img/8-points%20Decimate-in-time%20Algorithm.jpg">
>> </picture>
>> 
>> **Time Complexity :**\
>> Three layers of butterfly diagram, each layer with time complexity n \$\rightarrow\$ \$3\times n\$
>
> #### Radix-2 Cooley–Tukey Algorithm (Carl Friedrich Gauss' Algorithm) - \$\mathcal{O}(n\log{n})\$
>> $$
>> \begin{aligned}
>>   X_k
>>     &= \sum_{n\ even}^{N-1}\ \ x_n\ e^{-\frac{i2\pi}{N}nk} \qquad + \ \sum_{n\ odd}^{N-1}x_n\ e^{-\frac{i2\pi}{N}nk} \\
>>     &= \sum_{m=0}^{(N/2)-1}x_{2m}\ e^{-\frac{i2\pi}{N}(2m)k} \ + \sum_{m=0}^{(N/2)-1}x_{2m+1}\ e^{-\frac{i2\pi}{N}(2m+1)k} \\
>>     &= \sum_{m=0}^{(N/2)-1}x_{2m}\ e^{-\frac{i2\pi}{N/2}mk} \quad + e^{-\frac{i2\pi}{N}k}\sum_{m=0}^{(N/2)-1}x_{2m}\ e^{-\frac{i2\pi}{N/2}mk} \\
>>   X_{k+\frac{N}{2}} 
>>     &= \sum_{m=0}^{(N/2)-1}x_{2m}\ e^{-\frac{i2\pi}{N/2}m(k+\frac{N}{2})}\ + e^{-\frac{i2\pi}{N}k}\sum_{m=0}^{(N/2)-1}x_{2m+1}\ e^{-\frac{i2\pi}{N/2}m(k+\frac{N}{2})} \\
>>     &= \sum_{m=0}^{(N/2)-1}x_{2m}\ e^{-\frac{i2\pi}{N/2}mk} \qquad\ - e^{-\frac{i2\pi}{N}k} \sum_{m=0}^{(N/2)-1}x_{2m+1}\ e^{-\frac{i2\pi}{N/2}mk} \\
>> \end{aligned}
>> $$
>> 
>> Summary, 
>>
>> $$
>> \begin{aligned}
>>   X_k
>>     &= EvenPart_k + e^{-\frac{i2\pi}{N}k}OddPart_k \\
>>   X_{k+\frac{N}{2}}
>>     &= EvenPart_k - e^{-\frac{i2\pi}{N}k}OddPart_k
>> \end{aligned}
>> $$
>> 
>> **Time Complexity :**\
>> \$\log{n}\$ layers of butterfly diagram, each layer with time complexity n \$\rightarrow\$ \$(\log{n})\times n\$
>
> #### Radar's Algorithm - \$\mathcal{O}(n\log{n})\$
>> Rewrite Discrete Fourier Transform Algorithms
>> 
>> $$
>> \begin{align}
>>   X_0 & = \sum_{n=0}^{N-1}x_n \\ 
>>   X_{g^{-p}} & = x_0 + \sum_{q=0}^{N-2}x_{g^q}\ e^{-\frac{i2\pi}{N}g^{-(p-q)}} \qquad p = 0, \cdots, N-2 
>> \end{align}
>> $$
>> 
>> where n \$\in\$ { 1, 2, ..., N-1 } biject to q \$\in\$ { 0, 1, ..., N-2 }, k \$\in\$ { 1, 2, ..., N-1 } biject to p \$\in\$ { 0, 1, ..., N-2 }[^1].
>> 
>> $$ \sum_{q=0}^{N-2}x_{g^q}\ e^{-\frac{i2\pi}{N}}g^{-(p-q)} = x_{g^q\pmod{N}} \otimes e^{-\frac{i2\pi}{N}g^q\pmod{N}} $$
>> Assume
>> 
>> $$
>> \begin{array}{c}
>>   a = x_{g^q\pmod{N}} \quad\ \ \\
>>   b = e^{-\frac{i2\pi}{N}g^q\pmod{N}} 
>> \end{array} \qquad q = 0, \cdots, N-2
>> $$
>>
>> Completing the a-array and b-array to the smallest power of 2 greater than 2(N-1)-1[^2]:
>>  - a - padding zeros between first (index = 0) and second (index = 1) elements.
>>  - b - Repeat itself.
>> 
>> $$ \Large{\begin{align}
 \hat{a} \ast \hat{b} & = \sum_{m=0}^{M-1}\hat{a}\[m\]\hat{b}\[n-m\] \\
 & = a\[0\]\hat{b}\[n\] + \underbrace{0\cdot \hat{b}\[n-1\] + \cdots + 0\cdot \hat{b}\[n-M+N-1\]}_{\normalsize{M-N+1}} \\
 & \qquad\qquad + \underbrace{a\[1\]\hat{b}\[n-M+N-2\] + \cdots + a\[N-2\]\hat{b}\[n-M+1\]}_{\normalsize{N-2}} \\
 & = a\[0\]b\[n\] + a\[1\]b\[n-(M-N+2)_{mod\ N}\] + \cdots + a\[N-2\]b\[n-(M-1)_{mod\ N}\] \\
 & = \hat{a} \otimes \hat{b} = IDFT(\ DFT(\hat{a})\cdot DFT(\hat{b})\ )
 \end{align} }$$
>>
>> **Time Complexity :**\
>> Time complexity of DFT \$(\hat{a})\approx n\log{n}\$, and same as DFT \$(\hat{b})\$ and IDFT \$\rightarrow\$ \$3\times n\log{n}\$
> 
> #### (Hybrid Radix) Cooley–Tukey Algorithm - \$\mathcal{O}(n\log{n})\$
>> $$\Large{\underbrace{\begin{bmatrix}\ a_0, a_1, \cdots, a_{\normalsize N-1}\ \end{bmatrix}}_{\normalsize{N}} \Longrightarrow 
 \begin{bmatrix}
 a_{\normalsize 0,\ 0} & \cdots & a_{\normalsize 0,\ N2-1} \\
 \vdots   & \ddots & \vdots   \\
 a_{\normalsize N1-1,\ 0}& \cdots & a_{\normalsize N1-1,\ N2-1}
 \end{bmatrix}_{N1\times N2}}$$
>>
>> Rewrite Discrete Fourier Transform Algorithms
>> 
>> $$
>> \begin{align}
>>   X_{N_2k_1+k_2}
>>     &= \sum_{n_1=0}^{N_1-1}\sum_{n_2=0}^{N_2-1}x_{N_1n_2+n_1}\ e^{-\frac{i2\pi}{N_1N_2}(N_1n_2+n_1)(N_2k_1+k_2)} \\ 
>>     &= \underbrace{\sum_{n_1=0}^{N_1-1} (\ (\overbrace{e^{-\frac{i2\pi}{N_1N_2}n_1k_2}}^{\normalsize{Twiddle\ Factor}}) (\overbrace{\sum_{n_2=0}^{N_2-1}x_{N_1n_2+n_1}\ e^{-\frac{i2\pi}{N_2}n_2k_2}}^{\normalsize{Do\ N1\ times\ of\ DFT\ transformations\ with\ N2\ length}})\ )e^{-\frac{i2\pi}{N_1}n_1k_1}}_{\normalsize{Do\ N2\ times\ of\ DFT\ transformations\ with\ N1\ length}}
>> \end{align}
>> $$
>> 
>> **Time Complexity :**\
>> Time complexity of inner DFT is \$N_1\log{N_2}\ ;\$\
>> outer DFT is \$N_2\log{N_1}\$ \$\rightarrow\$ \$N_1\log{N_2}\times N_2\log{N_1} = N(\log{N_1} + \log{N_2})\$

[^1]: Properties of multiplicative group of integers modulo n.
[^2]: Since convolution theorem, convolution result is same as circular convolution result when origin array length is padding to more than 2(N-1)
