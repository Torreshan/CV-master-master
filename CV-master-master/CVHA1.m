%  Gruppennummer:M08
%  Gruppenmitglieder:Qu,Xingwei; Han,Fengze; Gao,shang; Chen,Zhong; Wu,
%  Yingyu

%% Hausaufgabe 1
%  Einlesen und Konvertieren von Bildern sowie Bestimmung von 
%  Merkmalen mittels Harris-Detektor. 

%  F? die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter ?er den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden k?nen.


%% Bild laden
  Image = imread('szene.jpg');
   IGray = rgb_to_gray(Image);
  
   
   
   
   
   

%% Harris-Merkmale berechnen
    tic;
     Merkmale = harris_detektor(IGray,6,0.05,0.1,10,[200,300]);
    %Merkmale = harris_detektor(IGray);

   
   %% the function is function [Merkmale] = harris_detektor(Image,segment_length,k,tau,min_dist,tile_size,N)
   %% and the default value of these varibles is segment_length=6; k=0.05;tau=0.1;min_dist=10;tile_size=[200,300];N=10;
   %% It's also ok, when you input only a few parameters.

    toc;
