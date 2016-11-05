!*************************************************************************
!******   CODAGE DE L'EQUATION D'ADVECTION LINEAIRE SCALAIRE   ***********
!*************************************************************************
!*************************************************************************
!***** Comment editer votre fichier .f90 (Quels editeurs de texte ?) *****
!*****     - sous Windows : Crimson Editor...                        *****
!*****     - sous UNIX/Linux : nedit, emacs...                       *****
!*************************************************************************
program advection
  implicit none

  !##### DECLARATION DES VARIABLES #####
  !***** Parametres *****       ! parameter : variables qui resteront constantes durant l'evolution du code
  integer, parameter :: nx1=101 ! integer : entier relatif ; nx1 : nombre de points du maillage en x

  !***** Generale *****
  integer :: i,j

  !***** Maillage  *****
  double precision :: delta_x1             ! double precision : reel / delta_x1 : pas d'espace
  double precision :: xg                   ! limite a gauche du domaine
  double precision :: xd                   ! limite a droite du domaine
  double precision :: long_x1              ! longueur du domaine
  double precision, dimension(1:nx1) :: x1 ! tableau de reels / discretisation du domaine

  !***** Discretisation en temps *****
  integer :: nt                             ! nombre de pas de temps
  double precision :: delta_t               ! pas de temps

  !***** Champ de quantite advectee *****
  double precision, dimension(1:nx1) :: w1_nplus1,w1_n  ! quantite advectee

  !***** Vitesse d'advection *****
  double precision :: adv                               ! vitesse d'advection

  !***** Autres parametres *****
  double precision :: theta                             ! parametre theta

  !##### DEBUT DU CODE #####
  !***** Initialisation *****
  !----- Domaine de calcul -----
  xg=-2.0D+00
  xd=+2.0D+00

  long_x1=xd-xg
  delta_x1=long_x1/(nx1-1)

  !+++++ Affichage a l'ecran du pas d'espace +++++
  !+++++ En FORTRAN, le '06' indique toujours une sortie sur l'ecran +++++
  write(06,*) 'Delta_x1 = ',delta_x1

  do i=1,nx1
     x1(i)=(i-1)*delta_x1+xg
  enddo




  !----- Champ de quantite advectee - Etat initial du champ de quantite advectee
  ! exp(y) = fonction exponentielle de y
  !----- Parametre du calcul -----
  adv= 0.5D+00
  delta_t=12D-02 !2D-02 !8D-02  !12D-02     ! respectivement Cas !CFL<1 !CFL=1 !Cas CFL>1
  nt=300

  
  do i=1,nx1
  W1_n(i)=exp(-0.5*(((x1(i)/0.4)*(x1(i)/0.4))))  !Initialisation de la gaussienne pour t=0
  enddo

  theta= 0D-02!(adv*(delta_t/delta_x1))*(adv*(delta_t/delta_x1)) ! valeur du theta pour l'amortissement numérique (stabilité)

  !***** Calcul de la solution *****
  !----- Ouverture du fichier de sortie -----
  !+++++ Attention : - en FORTRAN les choix pour unit sont implicites jusque '09' +++++
  !+++++               par exemple, le '06' correspond a la sortie sur ecran      +++++
  !+++++             - a une unite correspond un et un seul fichier de sortie     +++++
  !+++++               jusqu'a ce que celui-ci ait ete ferme par l'instruction    +++++
  !+++++               'close'                                                    +++++
  open(unit=10,file='resultats.dat')

  !----- En-tete de fichier pour Tecplot (a commenter pour Gnuplot) -----
  !write(10,*) 'TITLE="RESULTATS"'
  !write(10,*) 'VARIABLES="X","U"'
  !write(10,*) 'ZONE I=',nx1

  !+++++ Ecriture de la solution initiale dans le fichier de sortie +++++
  do i=1,nx1
     write(10,'(2(f8.4,1x))') x1(i),w1_n(i)
  enddo

  !+++++ Pour une sortie sur Gnuplot (a commenter pour Tecplot) +++++
  write(10,*)
  write(10,*)

  !----- Schema numerique -----
  do j=1,nt ! Avancement temporel (Boucle en temps)

     !+++++ Mise en place du schema de dicretisation en espace +++++
     do i=2,nx1 ! Attention au choix des bornes de la boucles : il ne faut pas deborder des dimensions du tableau
      ! w1_nplus1(i) =  w1_n(i) - ((adv*delta_t)/delta_x1)*(w1_n(i)-w1_n(i-1))  !schema 2 Upwind décentré
       w1_nplus1(i) = w1_n(i) - (((adv*delta_t)/(2*delta_x1))*(w1_n(i+1)-w1_n(i-1))) + ((theta/2)*(w1_n(i+1)-2*w1_n(i)+w1_n(i-1))) !schema 3 centré
      enddo


      do i=2,nx1
         w1_n(i)=w1_nplus1(i) ! remplacement de la valeur w1_n par la valeur w1_nplus1 afin de passer à l'itération suivante
      enddo
      

      w1_n(1) = w1_n(nx1) ! fermeture de la boucle permettant le calcul de la dernière valeurs (nx1) grâce à la premiere (qui pour elle est la                                  suivante)

     !+++++ Que faire des points que l'on n'a pas pu traiter = Conditions aux limites                 +++++
     !+++++ Dans la majorite des cas que vous aurez a etudier, on utilisera principalement            +++++
     !+++++ deux types de conditions limites :                                                        +++++
     !+++++ - les conditions de Dirichlet : w=0                                                       +++++
     !+++++ - les conditions de Von Neumann : dw/dn=0 , n etant la normale a la paroi                 +++++
     !+++++ Dans le projet qui vous est propose, on vous demande en plus des conditions periodiques : +++++
     !+++++ Astuce : recopier la sortie sur l'entree !                                                +++++

     !+++++ Ecriture des resultats dans fichier de sortie a chaque pas de temps +++++
     !+++++ Pour une sortie sur Tecplot (a commenter pour Gnuplot) +++++
     !write(10,*) 'ZONE I=',nx1
     !+++++ Pour une sortie sur Gnuplot (a commenter pour Tecplot) +++++
     write(10,*)
     write(10,*)

     !+++++ Ecriture +++++
     do i=1,nx1
        write(10,'(2(f8.4,1x))') x1(i),w1_n(i)
     enddo

  enddo

  !----- Fermeture du fichier de sortie -----
  close(10)

end program advection
!*****************************************************************************
!***** Pour plus d'information sur le langage FORTRAN 90/95 vous pouvez  *****
!***** telecharger la documentation sur le site de l'IDRIS a l'adresse : *****
!*****   http://www.idris.fr/data/cours/lang/fortran/choix_doc.html      *****
!*****   Choisir l'edition "Fortran : notions de base" (1er niveau)      *****
!*****************************************************************************
