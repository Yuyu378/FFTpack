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
> ### One-dimensional Discrete Fourier Transform
> #### Discrete Fourier Transform Algorithm - \$O(n^2)\$
> $$\Large{ X_k = \sum_{n=0}^{N-1}x_n\ e^{-\frac{i2\pi}{N}nk} \qquad for\ k = 0, \cdots, N-1 }$$

> ### One-dimensional Fast Fourier Transform
>> #### 8-points Decimate-in-time Algorithm - \$O(3n+5)\$
>> $$\Large{ W_n^0 = e^{-\frac{i2\pi}{N}0} = 1 }$$
>> 
>> ![8-points DIT](https://th.bing.com/th/id/R.dcef1a65d97e82adf10929a2e4eabfaf?rik=yrxhCGyT%2bGcuRA&riu=http%3a%2f%2fimg.blog.csdn.net%2f20170901201507165%3fwatermark%2f2%2ftext%2faHR0cDovL2Jsb2cuY3Nkbi5uZXQvZXpfeXd3%2ffont%2f5a6L5L2T%2ffontsize%2f400%2ffill%2fI0JBQkFCMA%3d%3d%2fdissolve%2f70%2fgravity%2fSouthEast&ehk=80flaSpbvyMnzt72FK6Ch6h6VNdG2MtpjumtaZQ36qE%3d&risl=&pid=ImgRaw&r=0)
>
>> #### Radix-2 Cooley–Tukey Algorithm - \$O(n\log{n})\$
>> $$\Large{\begin{align}
 X_k & = \ \sum_{n\ even}^{N-1}\ \ x_n\ e^{-\frac{i2\pi}{N}nk} \quad + \ \sum_{n\ odd}^{N-1}x_n\ e^{-\frac{i2\pi}{N}nk} \\
     & = \sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N}(2m)k} + \sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N}(2m+1)k} \\
     & = \underbrace{\sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N/2}mk}}_{\normalsize{DFT\ of\ even-indexed\ part\ of\ x_n}} + e^{-\frac{i2\pi}{N}k}\underbrace{\sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N/2}mk}}_{\normalsize{DFT\ of\ odd-indexed\ part\ of\ x_n}} \\
     \\
 X_{k+\frac{N}{2}} & = \sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N/2}m(k+\frac{N}{2})} + e^{-\frac{i2\pi}{N}k}\sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N/2}m(k+\frac{N}{2})} \\
     & = \underbrace{\sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N/2}mk}}_{\normalsize{DFT\ of\ even-indexed\ part\ of\ x_n}} - e^{-\frac{i2\pi}{N}k}\underbrace{\sum_{m=0}^{(N/2)-1}x_n\ e^{-\frac{i2\pi}{N/2}mk}}_{\normalsize{DFT\ of\ odd-indexed\ part\ of\ x_n}} \\
\end{align} }$$
>> 
>> Summary, 
>> 
>> $$\Large{\begin{align}
 X_k & = E\normalsize{venPart}\Large{_k + e^{-\frac{i2\pi}{N}k}O}\normalsize{ddPart}\Large{_k} \\
 X_{k+\frac{N}{2}} & = E\normalsize{venPart}\Large{_k - e^{-\frac{i2\pi}{N}k}O}\normalsize{ddPart}\Large{_k}
\end{align} }$$
>> 
