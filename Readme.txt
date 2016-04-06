Readme



ACF_figures_and_indexes_for_ACP:
Permet d'affichier une chronique avec son ACF et l'intervalle de confiance associ� et de les enregistrer
Permet de calculer les indices caract�ristiques de l'ACP sur l'ACF et de les enregistrer





Fourier_figure_and_indexes_for_ACP:
Permet d'affichier la FFT d'une chronique et son intervalle de confiance associ� et de l'enregistrer enregistrer
Permet de calculer les indices caract�ristiques  de l'ACP sur le spectre et de les enregistrer





ezfft
Realise la transformee de Fourier et affiche automatiquement le spectre de Fourier
En entree un vecteur et une chronique
En sortie [W,E]: les frequence et leur energie





ACF_Fourier_figures_and_indexes:
Code utilis� pour BVH
Calcul pour les donn�es brutes et les anoamlies
L'ACF et son intervalle de confiance
La transform�e de Fourier et son intervalle de confiance
Les indices caract�ristiques pour l'ACP

En sortie: une representation graphique avec
En colonne de droite; la chronique brute, son ACF et intervalle de confiance, son spectre de fourier et intervalle de confiance
En colonne de gauche; les anoamlies, son ACF avec intervalle de confiance, son spectre de fourier et intervalle de confiance
Les images sont automatiquement enregistr�es dans un repertoire dont on doit definir le nom

Deux fichiers textes dont on doit d�finir le nom et le repertoire
L'un donne les indices caract�ristiques de l'ACP pour l'ACF et le spectre de Fourier des donn�es brutes
L'autre donne les indices caract�ristiques de l'ACP pour l'ACF et le spectre de Fourier des anomalies
