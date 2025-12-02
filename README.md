## Het modelleren van tumoren met behulp van ODE's: ##


## Gebruikte wiskundige modellen en toelichting: ###
De groei van tumoren kan worden beschreven door middel van differentiaalvergelijkingen: de verandering in het volume 
($V$) van een tumor over de tijd ($t$) is een functie van het huidige volume van de tumor. $$
\frac{dV}{dt} = f(V)
$$
Hoewel het volume varieert over de tijd, $V(t)$, hangt de functie $f$ zelf niet expliciet af van de tijd. Het model is
tijd-invariant (d.w.z. het gedrag van de groei van een tumor hangt niet af van het moment waarop 
hij begonnen is te groeien). Afhankelijk van de keuze en aannames rond de functie $f$ kunnen verschillende modellen 
worden verkregen. Zoals onderstaande modellen.

### **lineair model:** 
$$
\frac{dV}{dt} = c
$$
Een lineaire groeifunctie beschrijft een hoeveelheid die met een constante snelheid in de loop van de tijd toeneemt 
of afneemt.
In dit model is de groeisnelheid $c$ constant en onafhankelijk van het volume $V$. 
  
**Toepassing bij tumorgroei:**  
In het geval van tumorgroei 
kun je zeggen dat het lineaire model tumorgroei met een constante groeisnelheid, onafhankelijk van de tumorgrootte 
beschrijft.  

**Voordelen/pluspunten:**
- Simpel en intuïtief model.  
- Makkelijk te analyseren en te kalibreren.  

**Beperkingen:**  
- Houd geen rekening met afvlakking van groei bij grotere tumoren.  
- Negeert beïnvloeding van de groeisnelheid door bijvoorbeeld voedingsstoffen of ruimte.


### **exponentiëel model:** 
$$
\frac{dV}{dt} = c \cdot V
$$
Een exponentiële groeifunctie beschrijft een groei waarbij de snelheid evenredig is met het huidige volume; hoe groter de tumor wordt, des te sneller hij groeit. 
De groeisnelheid $c$ neemt dus toe met toenemend volume $V$.

**Toepassing bij tumorgroei:**
In het geval bij tumorgroei kun je zeggen dat het exponentiële model een tumor beschrijft die groeit met 
een snelheid die evenredig is aan zijn volume.
Geschikt voor kleine tumoren in vroege stadia waarbij de remmende factoren minimaal zijn.

**Voordelen/pluspunten:**
- goed voor snelle groei
- analytisch eenvoudig oplosbaar

**Beperkingen:**  
- Groei is onbeperkt, dit is biologisch vaak onrealistisch.


### **Mendelsohn model**
$$
\frac{dV}{dt} = c \cdot V^d
$$
Een Mendelsohn model generaliseert exponentiële groei door een exponent $d$ toe te voegen. De waarde van deze exponent, 
 $d$, is vaak gerelateerd aan de geometrie van de groei en kan verschillende waarden aannemen, afhankelijk van het 
specifieke model of de fase van tumorgroei. Dit model stelt dat de groeisnelheid evenredig is met het volume tot een bepaalde macht, wat vaak resulteert in een groei die sneller is dan lineair, maar trager is dan exponentieel.

**Toepassing bij tumorgroei:**
kan de snelheid van tumorgroei realistisch aanpassen afhankelijk van de tumoromstandigheden.

**Voordelen/pluspunten:**
- Flexibel model dat groei kan vertragen $d$<1 of versnellen $d$>1

**Beperkingen:**
- Extra parameter maakt calibratie complexer.

### **exponentieel afvlakkend model** 
$$
\frac{dV}{dt} = c \cdot (V_{\max} - V)
$$
Bij het exponentieel afvlakkend model vertraagt de groei naarmate het volume dichter bij $V_{\max}$ komt. Bij dit model is de groeisnelheid evenredig met de resterende ruimte tot het maximum, de groei is in het begin het snelst en neemt continu af naarmate het volume het limiet nadert.

**Toepassing bij tumorgroei:**
Geschikt om verzadiging en limieten in tumorvolume te modelleren. 

**Voordelen/pluspunten:**
- Simpele manier om afvlakking in te bouwen

**Beperkingen:**  
- Vereist kennis van het maximale volume  $V_{\max}$

### **logistisch model**
$$
\frac{dV}{dt} = c \cdot V \cdot (V_{\max} - V)
$$
Beschrijft snelle initiële groei die afvlakt naarmate het volume $V_{\max}$ nadert. 

**Toepassing bij tumorgroei:** 
veel gebruikt om tumoren te modelleren die een draagcapaciteit hebben. Dit model beschrijft een S-curve: de tumor groeit initieel exponentieel, maar vertraagt lineair naarmate het volume de maximale hoeveelheid nadert door beperkte ruimte of voedingsstoffen.   

**Voordelen / pluspunten:**  
- Realistisch initieel exponentiële groei gevolgd door afvlakking.  

**Beperkingen:**  
- Vereist kennis van het maximale volume $V_{\max}$.  

## **Montroll model**
$$
\frac{dV}{dt} = c \cdot V \cdot (V_{\max}^d - V^d)
$$
Breidt het logistische model uit met een exponent $d$ voor flexibeler afvlakgedrag. Dit is een algemene vorm van het logistische en Gompertz-model, waarbij een extra parameter wordt toegevoegd om de vorm van de curve en het punt waarop de vertraging optreedt flexibeler te beschrijven.  

**Toepassing bij tumorgroei:** 
geschikt wanneer het afvlakken van groei niet lineair verloopt.  

**Voordelen / pluspunten:**  
- Realistischer afvlakgedrag dan standaard logistisch model.  

**Beperkingen:**  
- Extra parameter $d$ maakt modelcomplexiteit groter.


### **Allee-effect model**
$$
\frac{dV}{dt} = c \cdot (V - V_{\min}) \cdot (V_{\max} - V)
$$
Groei start pas boven een minimumvolume $V_{\min}$. Dit model houdt rekening met een drempelwaarde. 
 
**Toepassing bij tumorgroei:** 
nuttig voor tumoren die pas effectief groeien na een kritieke grootte. Als de tumorpopulatie te klein is, is de groeisnelheid negatief of nul (de tumor sterft uit), en groei treedt pas op zodra het volume boven deze minimale overlevingsgrens komt.  

**Voordelen / pluspunten:**  
- Beschrijft minimum effectieve tumoromvang.  

**Beperkingen:**  
- Vereist parameter $V_{\min}$.   


### **lineair gelimiteerd model**
$$
\frac{dV}{dt} = \frac{c \cdot V}{V + d}
$$
De groeisnelheid neemt af bij grotere volumes door de limiet $d$.  

**Toepassing bij tumorgroei:** 
geschikt voor tumoren met beperkte groeicapaciteit door ruimte of voedingsstoffen. Dit model beschrijft vaak tumoren waarvan de groeisnelheid begrensd wordt tot een constant maximum zodra de tumor een bepaalde grootte heeft bereikt. 

**Voordelen / pluspunten:**  
- Simpele afvlakking ingebouwd.  

**Beperkingen:**  
- Parameter $d$ nodig. 

### **oppervlakte-gelimiteerd model**
$$
\frac{dV}{dt} = \frac{c \cdot V}{(V + d)^{1/3}}
$$
De groei wordt beperkt door het oppervlak (diffusie van voeding).  

**Toepassing bij tumorgroei:** 
beschrijft fysiek gelimiteerde groei, realistisch voor grotere tumoren. Bij dit model groeien alleen de cellen aan de buitenrand van de tumor waardoor de groeisnelheid evenredig is met de oppervlakte in plaats van het totale volume.
 
**Voordelen / pluspunten:**  
- Houdt rekening met diffusielimieten.  

**Beperkingen:**  
- Parameter $d$ nodig.  

### **Von Bertalanffy model**
$$
\frac{dV}{dt} = c \cdot V^{2/3} - d \cdot V
$$
Metabolische verliezen verminderen de groeisnelheid naarmate het volume toeneemt.  

**Toepassing bij tumorgroei:** 
veel gebruikt in biologische groei-analyse. Dit model wordt ook gebruikt om de groei van individuele organismen te beschrijven.  

**Voordelen / pluspunten:**  
- Afvlakking door metabolisme meegenomen.  

**Beperkingen:**  
- Vereist parameters $c$ en $d$. 
- 
### **Gompertz model**
$$
\frac{dV}{dt} = c \cdot V \cdot \ln\left(\frac{V_{\max}}{V}\right)
$$
Initieel exponentiële groei die afvlakt richting $V_{\max}$. Vergelijkbaar met logistische groei, maar de afname van de specifieke groeisnelheid verloopt exponentieel.  

**Toepassing bij tumorgroei:** 
vaak gebruikt voor menselijke tumoren, realistisch verloop van initiële snelle groei tot verzadiging.  

**Voordelen / pluspunten:**  
- Realistisch voor veel biologische tumoren.  
- Beschrijft initieel snelle groei gevolgd door afvlakking.  

**Beperkingen:**  
- Vereist parameter $V_{\max}$. 


### 2. Biologische achtergrond: ###
Chinese Hamster V79 fibroblast tumor cellen worden veel gebruikt in onderzoek naar
DNA schade en DNA reparatie. De V79 p53 sequentie bevat twee punt mutaties in DNA bindende domeinen liggen. 
p53 is een tumorsuppressoreiwit dat betrokken is bij het beschermen van DNA en bestaat uit vier subunits.



### Bronnen: ###
1. Treatment Optimization for Tumor Growth by Ordinary Differential Equations
2. Euler's Method: https://tutorial.math.lamar.edu/classes/de/eulersmethod.aspx
3. Hamster cellen: https://pmc.ncbi.nlm.nih.gov/articles/PMC146528/
4. Groeimodellen: https://cdn.hackaday.io/files/8586367343424/From%20the%20Mendelsohn%20model%20to%20the%20Gompertz%20and%20logistic%20growth%20law.pdf
5. Groeimodellen: https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12450
6. Groeimodellen: https://set.kuleuven.be/voorkennis/blik-op-wiskunde/handboekB/explog/groeimodellen
7. Groeimodellen: https://cdn.hackaday.io/files/8586367343424/From%20the%20Mendelsohn%20model%20to%20the%20Gompertz%20and%20logistic%20growth%20law.pdf
8. Groeimodellen: https://www.sciencedirect.com/science/article/abs/pii/S0734975024000296

