## Het modelleren van tumoren met behulp van ODE's: ##


### 1. Gebruikte wiskundige modellen: ###

Lineare groei:
Bij een linear groeimodel verhoogd V met een constante hoeveelheid per tijdsstap.

Exponentieel toenemende groei:
Bij dit model is de groeisnelheid recht evenredig met het huidige volume; hoe groter de tumor wordt, des te sneller hij groeit.

Mendelsohn groei:
Dit model stelt dat de groeisnelheid evenredig is met het volume tot een bepaalde macht, wat vaak resulteert in een groei die sneller is dan lineair, maar trager dan exponentieel.

Exponentieel afvlakkende groei:
Bij dit model is de groeisnelheid evenredig met de resterende ruimte tot het maximum, de groei is in het begin het snelst en neemt continu af naarmate het volume het limiet nadert.

Logistische groei:
Dit model beschrijft een S-curve: de tumor groeit initieel exponentieel, maar vertraagt lineair naarmate het volume de maximale hoeveelheid nadert door beperkte ruimte of voedingsstoffen.

Montroll groei:
Dit is een algemene vorm van het logistische en Gompertz-model, waarbij een extra parameter wordt toegevoegd om de vorm van de curve en het punt waarop de vertraging optreedt flexibeler te beschrijven.

Allee effect groei:
Dit model houdt rekening met een drempelwaarde. Als de tumorpopulatie te klein is, is de groeisnelheid negatief of nul (de tumor sterft uit), en groei treedt pas op zodra het volume boven deze minimale overlevingsgrens komt.

Lineair gelimiteerde groei:
Dit model beschrijft vaak tumoren waarvan de groeisnelheid begrensd wordt tot een constant maximum zodra de tumor een bepaalde grootte heeft bereikt.

Oppervlakte-gelimiteerde groei:
Bij dit model groeien alleen de cellen aan de buitenrand van de tumor waardoor de groeisnelheid evenredig is met de oppervlakte in plaats van het totale volume.

Von Bertalanffy groei:
Dit model wordt ook gebruikt om de groei van individuele organismen te bescrijven.

Gompertz groei:
Vergelijkbaar met logistische groei, maar de afname van de specifieke groeisnelheid verloopt exponentieel.


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
