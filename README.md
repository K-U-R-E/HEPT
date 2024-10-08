# HEPT - Hybrid Engine Performance Tool

The code performs a simulation of a nitrous oxide - paraffin hybrid rocket motor. Some inputs from the user are needed, as seen below. 

Authors: 

- Vinay Williams and Mohammed Zweiri, 2020. (Kingston University London)

- Guilherme Tavares, 2019.


## ASSUMPTIONS

- Isentropic flow along the nozzle

- Temperature inside the combustion chamber does not change along the
  burn

- Combustion chamber is adiabatic

- Oxidiser tank empties adiabatically (and temperature does not change in the vapour-only phase)

- Available volume of combustion chamber does not change

- Combustion products form a mixture which behaves like an ideal gas

- Ideal burning (no erosion along the burning port)

- No O/F ratio change along the burn (strong assumption, but with N20 the performance variation is not so dependant on O/F)

## KNOW BEFORE RUNNING

- Some inputs should be provided by the user on the section below.

- Scipy and CoolProp needed.

- No prior computation about the ullage for the oxidiser tank is carried out by this version of the code. Carefully compute how much         oxidiser will be pumped inside the tank. Also, carefully compute the tank's diameter and length.

- The vent of the oxidiser tank is considered to be on the top of the tank.

- For plotting different variables, proceed to the bottom of the code and follow the examples.

- Step can be reduced for further precision, but must be accompanied by a reduction on "while" conditions.

- New versions of the code will be gradually uploaded to tackle the above-mentioned assumptions, and to improve user experience.

## REFERENCES (Some major ones)

- **For general internal ballistic's equations:** SUTTON, G.P.; BIBLARZ, Oscar. Rocket Propulsion Elements. 8th ed. John Wiley & Sons, Inc., Hoboken, New Jersey, 2010.

- **For self-pressurizing N20 Model:** WHITMORE, S.A.; CHANDLER, S.N. Engineering Model for Self-Pressurizing Saturated-N20-Propellant Feed Systems. Journal of Propulsion and Power, Vol. 26, No. 4, July-August 2010.

- **For N20-paraffin's regression law:** LESTRADE, J.-Y.; ANTHOINE, J.; LAVERGNE, G. Liquefying fuel regression rate modeling in hybrid propulsion. The French Aerospace Lab, Mauzac, France; The French Aerospace Lab, Toulouse, France, 2010.

- **For N20 properties:** CoolProp (available documentation in http://www.coolprop.org/index.html)

## LICENSE

MIT License

Copyright (c) 2020 Vinay Williams and Mohammed Zweiri

Copyright (c) 2019 Guilherme Castrignano Tavares

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

***Note:*** *The code uses lines written by Guilherme which was licensed using an MIT 
      license therefore the above license inclusion. Furthermore, in an attempt
      to maintain open use of this software at Kingston University and other educational institutions,
      Mohammed and Vinay have also opted to distribute the 2020 version under a MIT License as well.*
