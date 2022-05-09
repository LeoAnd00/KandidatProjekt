README
===============
**Den här GitHub-mappen är en del av ett kandidatarbete vid Chalmers tekniska högskola: MVEX01-22-03 | Dynamisk modellering av signalvägar som reglerar kol- och kvävemetabolismen i bagerijäst**

**Sammandrag**

Förståelse för hur signalerign fungerar inom och mellan celler är relevant för att förstå sig på flera sjukdomar. En central process inom alla organismer är hur deras signalnätverk reagerar baserat på ändringar av näringstillgångar som exempelvis tillgången på kol och kväve. Signalvägarna TORC1, SNF1 och cAMP-PKA i bagerijäst är centrala för reglering av kol- samt kvävemetabolismen och är därför av intresse. I detta projekt undersöks dessa signalvägar genom modellering där en dynamisk ODE-modell baserad på artikeln skriven av Jalihal m.fl (https://doi.org/10.1091/mbc.E20-02-0117) används. För att analysera modellen implementerades denna och figurerna i Jalihal m.fl. återskapades för att bekräfta implementeringen. Modellen består av ett flertal parametrar som parameteruppskattas med syftet att ge en bättre beskrivning av data. Detta görs genom användning av en kvasi-Newtonalgoritm samt beräkning av identifierbarhet och känslighet för parametrarna. 

Två parametervektorer erhölls som ger en bättre beskrivning av datan. Slutsatsen som drogs från identifierbarhets- och känlighetsanalysen var att majoriteten av parametrarna hade låg känslighet. Kvalitativ data används för att undersöka modellens giltighet för scenarion där ingen kvantitativ data finns tillgänglig. Detta gjordes för de parametervektorer som erhölls från parameteruppskattingen och Jalihal m.fl. (https://doi.org/10.1091/mbc.E20-02-0117). Resultatet från denna analys var att de flesta scenarion som analyserades överensstämde med kvalitativ data. Det fanns även scenarion där detta inte var fallet. För dessa fall föreslogs experiment som kan förbättra och utöka modellen. Trots de brister modellen uppvisar är den ändå användbara som en grund för vidare utveckling. Åtgärder som utvidgning av modellen kombinerat med nya parameteruppskattningar kan på sikt leda till en bra beskrivande modell över signalvägarna i bagerijäst.

## GitHub-mappens struktur

- **Code** innehåller alla koder för att utföra bland annat återskapning av figurer från artikeln skriven av Jalihal m.fl. (https://doi.org/10.1091/mbc.E20-02-0117) samt parameteruppskattning.

- **Data** innehåller den experimentella datan som använts i detta projekt.

- **Intermediate** innehåller intermediära steg under parameteruppskattningen.

- **Results** innehåller resultaten från bland annat parameteruppskattningen och alla återskapade figurer. 

## Nödvändig programvara

Detta projekt genomfördes med programmeringsspråket **Julia**, version 1.7, inklusive vissa paket som kan ses i Project.toml filen.

# Skriven av:
- Leo Andrekson
- Emma Johansson
- Rakel Hellberg
- Jesper Olsson
 
