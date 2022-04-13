# -*- coding: utf-8 -*-
import math
import numpy as np

class Transformacje:
    
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            półosie:
            a - duża- promień równikowy
            b - mała - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.e = math.sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.e2 = (2 * self.flat - self.flat ** 2) # eccentricity**2

    def hirvonen(self,X, Y, Z):
        '''
        Użytkownik podaje następujące argumanety (współrzędne geocentryczne):
        X : wartość współrzędnej x danego punktu
        Y : wartość współrzędnej y danego punktu
        Z : wartość współrzędnej z danego punktu 
    
        Po użyciu funkcji użytkownik otrzymuje następujące wyniki (współrzędne geodezyjne):
        φ (fi)     : szerokość geograficzna punktu
        λ (lambda) : długość geograficzna punktu 
        H : wysokość punktu
        
       
        '''
       
        r = np.sqrt(X**2 + Y**2)
        fi_n = math.atan(Z / (r*(1-self.e2)))
        eps = 0.000001/3600*np.pi/180
        fi = fi_n * 2
        while math.sqrt((fi_n - fi)**2) > eps:
            fi = fi_n
            N = self.a/math.sqrt(1 - self.e2*math.sin(fi_n)**2)
            h = r/math.cos(fi_n) -N
            fi_n = math.atan(Z/(r*(1-self.e2*(N/(N+h)))))
    
        N = self.a / math.sqrt(1 - self.e2 * math.sin(fi_n) ** 2)
        h = r / math.cos(fi_n) - N
        lam = math.atan(Y/X)
    
    
        return fi_n, lam, h

  

    



    def uklad2k(self, fi, lam):
        '''
        Użytkownik podaje następujące argumenty:
            
        1) φ (fi)     :  szerokość geograficzna 
        2) λ (lambda) : długość geograficzna
    
        
        Po zostosowaniu funkcji uklad2k użytkownik otrzymuje następujące współrzędne :
            
        x00 : jest to współrzędna x przetransformowana do układu płaskiego 2000
        y00 : jest to współrzędna y przetransformowana do układu płaskiego 2000 
    
    
        '''
       
        m_0 = 0.999923
        N = self.a/(math.sqrt(1-self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = self.e2 * np.cos(lam)**2
        lam = math.degrees(lam)
        
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24
            
        lam = math.radians(lam)
     
        lam_0 = math.radians(lam_0)
        l = lam - lam_0
        
        A_0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A_2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A_4 = 15/256 * (self.e2**2 + (3*(self.e2**3))/4)
        A_6 = (35*(self.e2**3))/3072
        
        
        sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (s*1000000) + 500000, 3)   
        
        return x00, y00 

    def filamxyz(self, fi, lam,H):
        '''
        Użytkownik jako argumenty podaje następujące wartości :
            
        1) φ (fi)     :  szerokość geograficzna 
        2) λ (lambda) : długość geograficzna

        
        W wyniku działania funkcji użytkownik otrzymuje: 
            
        współrzędne geocentryczne : X, Y, Z
        
        '''
        
        N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
        
    
        X = round((N + H)*np.cos(fi)*np.cos(lam), 3)

        Y = round((N + H)*np.cos(fi)*np.sin(lam), 3)
  
        Z = round((N * (1 - self.e2) + H) * np.sin(fi), 3)

        
        return X, Y, Z

    def neu (self,x1, y1, z1, x, y, z):
        '''
        Użytkownik podaje następujące argumenty :
        
        Są to współrzędne geocentryczne punktu A oraz punktu B
        
       współrzędne x,y,z ----> punkt A
    
       współrzędne x,y,z ----> punkt B
        
        
        W wyniku działania funkcji użytkownik otrzymuje: 
            
        wartości wspolrzednych w trojwymiarowym układzie odniesienia : n,e,u
        
        '''
        fi, lam, h = self.XYZ2filamh(x, y, z)
        
        
        d_x = x - x1
        d_y = y - y1
        d_z = z - z1
        
        
        R = np.matrix([((-math.sin(fi) * math.cos(lam)), (-math.sin(fi) * math.sin(lam)), (math.cos(fi))),
                    ((-math.sin(lam)), (math.cos(lam)), (0)),
                    ((math.cos(fi) * math.cos(lam)), (math.cos(fi) * math.sin(lam)), (math.sin(fi)))])


        d = np.matrix([d_x, d_y, d_z])
        d= d.T
        neu = R*d
        N = neu[0,0]
        E = neu[1,0]
        U = neu[2,0]
        return N, E, U
    
    
    def uklad92(self, fi, lam ):
        '''
        Użytkownik podaje następuje następujące wartości:
            
        1) φ (fi)     :  szerokość geograficzna 
        2) λ (lambda) : długość geograficzna
    
        Po zostosowaniu funkcji uklad92 użytkownik otrzymuje następujące współrzędne :
            
        x92 : jest to współrzędna x przetransformowana do układu płaskiego 1992
        y92 : jest to współrzędna y przetransformowana do układu płaskiego 1992
    
        '''
    
        
        m_0 = 0.9993
        N = self.a/(math.sqrt(1-self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = self.e2 * np.cos(lam)**2
        lam_0 = math.radians(19) #poczatek ukladu w punkcie przeciecia poludnika L0 = 19st z obrazem równika 
        l = lam - lam_0
        
        A_0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A_2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A_4 = 15/256 * (self.e2**2 + (3*(self.e2**3))/4)
        A_6 = (35*(self.e2**3))/3072
        
        sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
        
        return x92, y92 
    
    def 2wymiar(self, xa, xb ,ya, yb):
        '''
      Użytkownik podaje następuje następujące wartości:
          
     Są to współrzędne punktu A oraz punktu B
          
        xa,ya ------> punkt A
        xb,yb ------> punkt B
    
       Po zostosowaniu funkcji użytkownik otrzymuje:
           
        odl czyli odegłość 2D
        '''
        
        odl = math.sqrt((xb - xa)**2 + (yb - ya)**2)
        return odl
    def 3wymiar(self, xa, xb, ya, yb, za, zb):
        odl = math.sqrt((xb - xa)**2 + (yb - ya)**2 +(zb-za)**2)
        return odl
        '''
      Użytkownik podaje następuje następujące wartości:
     
     Są to współrzędne punktu A oraz punktu B
     
        xa,ya,za ------> punkt A
        ya,yb,zb ------> punkt B
  

     Po zostosowaniu funkcji użytkownik otrzymuje:
      
     odl czyli odegłość 3D
        '''

    def AzEle(self,x1, y1, z1, x, y, z):
        '''
        Użytkownik podaje następuje następujące wartości:
            
    
          współrzędne x,y,z    ------> punkt A
            
          współrzędne x1,y1,z1 ------> punkt B
      
        Po zostosowaniu funkcji użytkownik otrzymuje:
            
            1) wartość wartość kąta azymutu (Az)
            2) wartość kąta elewacji (elewacja)
      
        '''
        
        N, E, U = Transformacje.NEU(x1, y1, z1, x, y, z)
        hz = math.sqrt(E**2 + N**2)
        Az= math.atan2(E,N)
        elewacja = math.atan2(U, hz)
        if Az <0:
            Az += 2* math.pi
        
        return Az, elewacja
        
















