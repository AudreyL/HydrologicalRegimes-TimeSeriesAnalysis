# HydrologicalRegimes-TimeSeriesAnalysis

This respository present the data and part of the code used in:
Gaucherel C., Frelat R. Lustig A., Rouy B., Chery Y. and Hubert P. (2016) Time–frequency analysis to profile hydrological regimes: application to Haiti, Hydrological Sciences Journal, vol 61., No. 2, 274-288.


Comprehensive characterization of flow rates is a prerequisite to a proper understanding andwater management of a given hydrological region. Several studies question the soundness of stationarity in time series and suggest the need for a quantification of the events and non-stationary features in flow rate time series. In this study, we combine statistical and time–frequency (TF) analyses to characterize and classify the flow rates of an understudied region, namely Haiti. Wavelet transforms and cyclostationarity analyses were combined with principal component analysis and hierarchical clustering to identify three groups of hydrological regimes in the country, suggesting similar management: (1) relatively stable flow rates with TF behaviour; (2) periodic and strongly seasonal flow rates; and (3) unstable flow rates. We argue that the TF methodology can yield additional information in regard to flow events and multiscale behaviour, even for short records. Flow rate characterization would benefit from the exhaustive approach described here. 

The data/code presented here include:
- ACF_and_indexes_for_ACP.m : A matlab code to display the autocorrelation Function of a streamflow and its confidence interval. It also calculate a set of 34 commonly used hydrological indices to characterize streamflows. 
- Fourier_and_indexes_for_ACP.m : A matalab code to display the Fourier Transform of a streamflow and its confidence interval. It also calculate a set of 34 commonly used hydrological indices to characterize streamflows.
- ezfft.m : Calculate the Fourier Transform
- ACF_Fourier_and_indexes : The matlab code used in the article to generate the Hydrological Profile (ACF, FOurier analysis and quantitative measures) of 24 streams in Haiti
- SSA.m : The Matlab code used to calculate the Singular Spectrum Analysis (SSA)
- Data : The flow rate statistics used in this study were provided by LGL (Lalonde, Girouard & Letendre Canadian consultants) S.A. freelance experts and the Haitian public services. The selected data set contained 26 daily time series pertaining to 17 watersheds for the 1919-1943 time period.
