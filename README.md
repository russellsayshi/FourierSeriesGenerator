# FourierSeriesGenerator
Generator for fourier series using my math parse repo. May extend it in a bit to solve the heat equation.

# How to use
Fire up main.cpp and it'll ask you about the series you want to generate, then spit out all the data ready for desmos. It prints in out the doubles fixed, which you can change in the code if you'd like. Here's some sample input/output:
```
Enter the function to create the fourier series (e.g. 2*sin(x)): x^2
What bound do you want on the function (Fourier series will be generated from -b to b): 5
How many terms do you want in your fourier series? 25
How many terms do you want in your Riemann sum? 100000
Cosine func: (x^2)*cos(n*pi*x/b)
Sine func: (x^2)*sin(n*pi*x/b)
8.3330833400100168262 + (-10.1316183775956627500)*cos(x*0.6283185307179586232) + (-0.0000000157074295541)*sin(x*0.6283185307179586232) + (2.5325296044327680889)*cos(x*1.2566370614359172464) + (0.0000000314147223615)*sin(x*1.2566370614359172464) + (-1.1252909427332933312)*cos(x*1.8849555921538758696) + (-0.0000000471219408006)*sin(x*1.8849555921538758696) + (0.6327574111313503114)*cos(x*2.5132741228718344928) + (0.0000000628293125862)*sin(x*2.5132741228718344928) + (-0.4047847479338040344)*cos(x*3.1415926535897931160) + (-0.0000000785366837239)*sin(x*3.1415926535897931160) + (0.2809477457056884231)*cos(x*3.7699111843077517392) + (0.0000000942440208547)*sin(x*3.7699111843077517392) + (-0.2062779391692419750)*cos(x*4.3982297150257103624) + (-0.0000001099513099701)*sin(x*4.3982297150257103624) + (0.1578143628171481683)*cos(x*5.0265482457436689856) + (0.0000001256586616412)*sin(x*5.0265482457436689856) + (-0.1245878944275353034)*cos(x*5.6548667764616276088) + (-0.0000001413660023326)*sin(x*5.6548667764616276088) + (0.1008211970417420000)*cos(x*6.2831853071795862320) + (0.0000001570733167589)*sin(x*6.2831853071795862320) + (-0.0832365288185412389)*cos(x*6.9115038378975439670) + (-0.0000001727806608318)*sin(x*6.9115038378975439670) + (0.0698619465099921677)*cos(x*7.5398223686155034784) + (0.0000001884879712284)*sin(x*7.5398223686155034784) + (-0.0594533765392665642)*cos(x*8.1681408993334621016)
```
And here is the lovely graph of that series:
![25 term fourier series of x^2](https://i.imgur.com/S0L8lZE.png)

Here is the same graph, expanded to 150 terms:
![100 term fourier series of x^2](https://i.imgur.com/yYf5QN3.png)

Almost indistinguishable from the real thing, isn't it! Ignore the little dips on the edges of the photo, those should be there. I just snipped the picture a little past -5 and 5, which were my bounds, and the series won't converge to the function outside of the bounds given.
