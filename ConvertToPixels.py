class ConvertToPixels:
    ChipSize = np.array([30., 15.]) # ALPIDE is a 30 mm x 15 mm chip
    PixelSize = np.array([0.02924,0.02688]) # pixelsize on ALPIDE
    # board contains 9 x 12 chips. The 12 chips are positioned in double layer 6-6 chips, with overlap overY \approx 0.1mm
    padX = 0.02912 # X padding: 29.12 micron
    Xgap = 0.1 # gap between ALPIDEs in X direction

    padY = 0.02944 # Y padding: 29.44 micron
    Ygap = 27.4 # distance between stripes in Y direction: 27.4 mm
    overY = 0.092 # overlap between F/B: 0.092 mm
    ElectrY = 1.208 # electronics part on the chip 1.208 mm
    Yshift = 0.85 # in Y direction the first chip starts at 0.85 mm (positive half)
    Ydif = 0.1 # for the negative half plane the first chip starts at 0.95 mm = 0.85 + Ydif
    
    ChipX = np.array([-ChipSize[0]/2+padX, ChipSize[0]/2-padX])
    ChipFY = np.array([Yshift + padY,Yshift + ChipSize[1] -ElectrY])
    YSize = (ChipFY[1] - ChipFY[0])/2. # 'radius' of the sensitive area (mm)
    ChipMidYF = (ChipFY[0]+ChipFY[1])/2
    ChipBY = np.array([Yshift + ChipSize[1] -ElectrY - overY + padY,Yshift + 2*ChipSize[1] -2*ElectrY - overY])
    ChipMidYB = (ChipBY[0]+ChipBY[1])/2.

    Xmin = 0
    Xmax = 1023
    Ymin = 0
    Ymax = 511
    
    Xmed = np.zeros(9) # middle of chip X positions, 9 chips 
    YmedF = np.zeros(6)
    YmedB = np.zeros(6)

    # cluster parameters (circleX, circleY Short_t circleX/Y[70])
    circleX = [0,1,0,-1,0,1,-1,-1,1,0,-2,0,2,1,-2,-1,2,-1,-2,1,2,-2,-2,2,2,0,-3,0,3,-1,-3,1,3,1,-3,-1,3,0,-4,0,4,2,-3,-2,3,-8,-2,-3,2,4,-1,-4,1,4,1,7,-1,3,3,-3,-3,4,2,-4,-2,4,-2,2,5,0];
    circleY = [0,0,-1,0,1,-1,-1,1,1,-2,0,2,0,-2,-1,2,1,-2,1,2,-1,-2,2,2,-2,-3,0,3,0,-3,1,3,-1,-3,-1,3,1,-4,0,4,0,-3,-2,3,2,15,-3,2,3,-1,-4,1,4,1,-4,-18,4,3,-3,-3,3,2,-4,-2,4,-2,-4,4,0,5];
    binPosLUT = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096]

    def __init__(self):
        self.Xmed = np.linspace(-4.*(self.ChipSize[0]+self.Xgap),4.*(self.ChipSize[0]+self.Xgap),9)
        self.YmedF = np.array((-83.15+self.padY-83.15+self.ChipSize[1]-self.ElectrY)/2+np.arange(6)*27.4)
        self.YmedF[3:] += 2*self.Yshift+self.Ydif
        self.YmedB = self.YmedF + 27.4/2
        self.PosZ = np.concatenate([[225.219,225.219+52.4],333.369+5.5*np.arange(0.,48.,1.)])  # Z position of Layers
    
    def get_Pos(self,X,Y,Edep):
        res = np.empty((2,4),dtype=int); res.fill(-1)
        col = np.argmin(np.abs(X-self.Xmed)) # selects chip (column)
        pX = (X+self.ChipX[1]-self.Xmed[col])/self.PixelSize[0]
        iX = np.floor(pX) # should be between 0 and 1023 for valid hit

        df = pd.DataFrame(columns=['col','row','X','Y','edep'],dtype=float)
        if(iX >= self.Xmin and iX <= self.Xmax):
            row = np.argmin(np.abs(Y-self.YmedF)) # selects chip (column)
            pY = (self.YSize+(Y-self.YmedF[row]))/self.PixelSize[1]  
            iY = np.floor(pY)  # should be between 0 and 511 for valid hit
            if(iY >= self.Ymin and iY <= self.Ymax):
                df = df.append({'col': col, 'row': 2*row, 'X': pX, 'Y': pY, 'edep': Edep}, ignore_index=True)
            row = np.argmin(np.abs(Y-self.YmedB)) # selects chip (column)
            pY = (self.YSize+(Y-self.YmedB[row]))/self.PixelSize[1]
            iY = np.floor(pY)  # should be between 0 and 511 for valid hit
            if(iY >= self.Ymin and iY <= self.Ymax):
                df = df.append({'col': col, 'row': 2*row+1, 'X': pX, 'Y': pY, 'edep': Edep}, ignore_index=True)
        return df
    
    def get_Pos_multi(self,ps):  # receives a list of X,Y pairs
#        return np.vstack([y for y in [self.get_Pos(p[0],p[1]) for p in ps] if 0 not in y.shape])
        df = pd.DataFrame(columns=['col','row','X','Y','edep'],dtype=float)
        for p in ps:
            df = df.append(self.get_Pos(p[0],p[1],p[2]))
        return df

    def get_Coords(self,p):
        X = p[2]*self.PixelSize[0] + self.Xmed[p[0]] - self.ChipX[1]
        j = p[1] // 2; back = p[1] % 2
        Ymed = self.YmedF[j]
        if(back == 1):
            Ymed = YmedB[j]
        if(j < 3):
            Y = Ymed + self.YSize - p[3]*self.PixelSize[1]
        else:
            Y = Ymed - self.YSize + p[3]*self.PixelSize[1]
        return [X,Y]
    
    def set_CS_multi(self,ps): # ps is [pX,pY,Edep]
        df = pd.DataFrame(columns=['iX','iY'],dtype=int)
        for p in ps:
            CS = np.floor(4.2267 * (p[2]*40.)**0.65 + 0.5)
            if(CS < 2): # cluster sizes less than 1 are under threshold
                continue
            if(CS < 27):  # use library for shapes
                id = np.random.randint(CSindex['start'][CS-1],CSindex['start'][CS])
                x_mean = CSconfigs['x_mean'][id]
                y_mean = CSconfigs['y_mean'][id]
                CSarray = CSconfigs['hit_array'][id]
                for i in range(10):              
                    for j in range(10):
                        if(CSarray[i] & binPosLUT[j]):
                            outX = np.floor(p[0] + int(np.floor(i-x_mean+0.5)))
                            outY = np.floor(p[1] + int(np.floor(j-y_mean+0.5)))
                            if(outX >= self.Xmin and outX <= self.Xmax and outY >= self.Ymin and outY <= self.Ymax):
                                df=df.append({'iX': outX, 'iY': outY}, ignore_index=True)
                continue
            CS = int(min(CS,70))  # max cluster size is 70
            for i in range(CS):
                outX = np.floor(p[0] + self.circleX[i]+0.5)
                outY = np.floor(p[1] + self.circleY[i]+0.5)
                if(outX >= self.Xmin and outX <= self.Xmax and outY >= self.Ymin and outY <= self.Ymax):
                    df=df.append({'iX': outX, 'iY': outY}, ignore_index=True)
        return df
    
    def get_Layer(self,Zs):
        return [np.sum(z > self.PosZ)-1 for z in Zs]
    
    def get_Layer_Pos(self):  # return the Z positions of Layers
        return self.PosZ
