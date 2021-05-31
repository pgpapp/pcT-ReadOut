// D:\root_v6.24.00\bin\root -b -q "C:/Users/Gábor/OneDrive - elte.hu/Documents/Research/pcT/ReadOut/get_database_root.cpp (\"DataFile\",\"CLusterDir\")"
// MicroSoft: D:\root_v6.24.00\bin\root -b -q 'C:/Users/Gábor/OneDrive - elte.hu/Documents/Research/pcT/ReadOut/get_database_root.cpp (""Data/MC/head_phantom_10000Primaries_230MeV_1617471898.root"",""Data/Cluster/"")'

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <array>
#include <random>

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iterator>

using namespace std;

class ConvertToPixels {
    private:
    std::array<int,2> ChipSize = {{30, 15}};  // ALPIDE is a 30 mm x 15 mm chip
    std::array<float,2> PixelSize = {{0.02924f,0.02688f}}; // pixelsize on ALPIDE
    // board contains 9 x 12 chips. The 12 chips are positioned in double layer 6-6 chips, with overlap overY \approx 0.1mm
    float padX = 0.02912f; // X padding: 29.12 micron
    float Xgap = 0.1f; // gap between ALPIDEs in X direction

    float padY = 0.02944f; // Y padding: 29.44 micron on the opposite to eletronics end
    float Ygap = 27.4f; // distance between stripes in Y direction: 27.4 mm
    float overY = 0.092f; // overlap between F/B: 0.092 mm on both ends
    float ElectrY = 1.208f; // electronics part on the chip 1.208 mm (only one side, no padding here)
    float Yshift = 0.85f; // in Y direction the first chip starts at 0.85 mm (positive half)
    float Ydif = 0.1f; // for the negative half plane the first chip starts at 0.95 mm = 0.85 + Ydif
    
    std::array<float,2> ChipX = {{-ChipSize[0]/2.f+padX, ChipSize[0]/2.f-padX}};
    std::array<float,2> ChipFY = {{Yshift + padY,Yshift + ChipSize[1] -ElectrY}};
    float YSize = (ChipFY[1] - ChipFY[0])/2.f; // 'radius' of the sensitive area (mm)
    float ChipMidYF = (ChipFY[0]+ChipFY[1])/2.f;
    std::array<float,2> ChipBY = {{Yshift + ChipSize[1] -ElectrY - overY + padY,Yshift + 2.f*ChipSize[1] -2.f*ElectrY - overY}};
    float ChipMidYB = (ChipBY[0]+ChipBY[1])/2.f;

    int Xmin = 0;
    int Xmax = 1023;
    int Ymin = 0;
    int Ymax = 511;

    std::array<float,9> XLow = {{0.f}};// middle of chip X positions, 9 chips 
    std::array<float,6> YLowF = {{0.f}};//  middle of "front" Y positions: 6 chips
    std::array<float,6> YLowB = {{0.f}};// middle of "back" Y positions: 6 chips

    // cluster parameters (circleX, circleY Short_t circleX/Y[70])
    // See Helge's https://github.com/HelgeEgil/DigitalTrackingCalorimeterToolkit repo and Thesis
    std::array<int,70> circleX = {{0,1,0,-1,0,1,-1,-1,1,0,-2,0,2,1,-2,-1,2,-1,-2,1,2,-2,-2,2,2,0,-3,0,3,-1,-3,1,3,1,-3,-1,3,0,-4,0,4,2,-3,-2,3,-8,-2,-3,2,4,-1,-4,1,4,1,7,-1,3,3,-3,-3,4,2,-4,-2,4,-2,2,5,0}};
    std::array<int,70> circleY = {{0,0,-1,0,1,-1,-1,1,1,-2,0,2,0,-2,-1,2,1,-2,1,2,-1,-2,2,2,-2,-3,0,3,0,-3,1,3,-1,-3,-1,3,1,-4,0,4,0,-3,-2,3,2,15,-3,2,3,-1,-4,1,4,1,-4,-18,4,3,-3,-3,3,2,-4,-2,4,-2,-4,4,0,5}};
    std::array<int,13> binPosLUT = {{1,2,4,8,16,32,64,128,256,512,1024,2048,4096}};

    std::array<float,50> PosZ = {{0.f}};
    float tol = 0.001f; // tolerance of positioning layer to posZ

    struct CS_struct {
        float x_mean;
        float y_mean;
        int size;
        std::array<int,10> bits;
    };

    std::vector<CS_struct> CSconfigs;
    std::vector<int> CSindex; // ranges for different clustersizes CSindex[i] is the last position of clustersize i in CSConfigs

    public:
        std::random_device rd{};
//       std::mt19937 generator(42);
        std::mt19937 generator{rd()};

        ConvertToPixels(std::string& CDIR) {
            float delta = ChipSize[0]+Xgap;
            for(size_t i=0; i<XLow.size(); i++) {XLow[i] = -4.f*(ChipSize[0]+Xgap) + ChipX[0] + i*delta;}
            for(int i=0; i<6; i++) {YLowF[i] = i*27.4f+(-83.15f+padY-83.15f+ChipSize[1]-ElectrY)/2.f-YSize;}
            for(int i=3; i<6; i++) {YLowF[i] += 2.f*Yshift+Ydif;}
            for(int i=0; i<6; i++) {YLowB[i] = YLowF[i] + 27.4f/2.f;}

            for(int i=0; i<2;i++) {PosZ[i] = 225.219f+i*52.4f;}
            for(int i=2; i<50;i++) {PosZ[i] = 333.369f+5.5f*(i-2);}  // Z position of Layers

			TFile        * CDB_fCluster; 
			CDB_fCluster = new TFile("C:/Users/Gábor/OneDrive - elte.hu/Documents/Research/pcT/ReadOut/Data/Cluster/database_final_reduced.root", "READ");
			TTreeReader reader("database", CDB_fCluster);
			TTreeReaderValue<float> x_mean(reader, "x_mean");
			TTreeReaderValue<float> y_mean(reader, "y_mean");
			TTreeReaderValue<Int_t> size(reader, "size");
			TTreeReaderArray<uint8_t> hits(reader, "hit_array");

			CS_struct t;
			while(reader.Next()) {
				t.x_mean = *x_mean;
				t.y_mean = *y_mean;
				t.size = *size;
				std::copy(hits.begin(),hits.end(),t.bits.begin());
				CSconfigs.push_back(t);
			}

            std::ifstream fin(CDIR+"\\sortIndex.csv");
            if(!fin) {cout << "Cannot open input file " << CDIR << "\\sortIndex.csv\n"; return;}
            int dummy;
            while (!fin.eof()) {
                fin >> dummy; fin >> dummy;
                CSindex.push_back(dummy);
            }
            fin.close();
        }
 

        std::vector<std::array<float,4>> get_Pos(float X,float Y,float edep) {
            std::vector<std::array<float,4>> res;
            int col = std::lower_bound(XLow.begin(), XLow.end(), X)-XLow.begin()-1; // select chip column
            if(col < 0 || col  >= (int)XLow.size()) return(res);
            float pX = (X-XLow[col])/PixelSize[0];
//            cout << "getPos: X=" << X << ", Y=" << Y << ", col=" << col << ", pX=" << pX << ", XLow[col] = " << XLow[col] <<endl;
            if(pX >= Xmin && pX <= Xmax) {
                int row = std::lower_bound(YLowF.begin(), YLowF.end(), Y)-YLowF.begin()-1; // select chip row
                if(row >= 0 && row < (int)YLowF.size()) {
                    float pY = (Y-YLowF[row])/PixelSize[1];
//                    cout << "\tgetPos1: row=" << row << ", pY=" << pY <<endl;
                    if(pY >= Ymin && pY <= Ymax) {
                        res.push_back({{1.f*col,2.f*row,pX,pY}});
                    }
                }
                row = std::lower_bound(YLowB.begin(), YLowB.end(), Y)-YLowB.begin()-1; // select chip row
                if(row >= 0 && row < (int)YLowF.size()) {
                    float pY = (Y-YLowB[row])/PixelSize[1];
//                    cout << "\tgetPos2: row=" << row << ", pY=" << pY <<endl;
                    if(pY >= Ymin && pY <= Ymax) {
                        res.push_back({{1.f*col,2.f*row+1.f,pX,pY}});
                    }
                }
            }
            return(res);
        }

        std::vector<pair<int,int>> get_CS_Pos(float pX, float pY, float edep) {  // position of activated pixels
            std::vector<pair<int,int>> res;
            int CS = floor(4.2267f*pow(edep*40.f,0.65)+0.5f);
            if(CS < 2) return(res);

            if(CS < 26) { // use library for shapes. THIS WAS 27 in Helges's code, however CSindex[26] == CSindex[25] ?!
                std::uniform_int_distribution<> dist(CSindex[CS],CSindex[CS+1]);
                int id = dist(generator);
                float x_mean = CSconfigs[id].x_mean;
                float y_mean = CSconfigs[id].y_mean;

                 for(size_t i=0; i<10; i++) {
                    for(size_t j=0; j<10; j++) {
                        if( CSconfigs[id].bits[i] & binPosLUT[j]) {
                            int iX = floor(pX + i - x_mean + 0.5);
                            int iY = floor(pY + j - y_mean + 0.5);
                            if(iX >= Xmin && iX <= Xmax && iY >= Ymin && iY <= Xmax) {
                                res.push_back(std::make_pair(iX,iY));
                            }
                       }
                    }
                }
                return(res);
            }

            CS = min(CS,70);
            for (int i=0; i< CS; i++) {
                int iX = floor(pX + circleX[i] + 0.5f);
                int iY = floor(pY + circleY[i] + 0.5f);
                if(iX >= Xmin && iX <= Xmax && iY >= Ymin && iY <= Xmax) {
                    res.push_back(std::make_pair(iX,iY));
                }

            } 

			return(res);
        }

#define SENSOR_THICKNESS 0.025f		
		int get_Layer(float z) {  // get the Layer id
			const auto& p = std::upper_bound(PosZ.begin(),PosZ.end(),z-tol-SENSOR_THICKNESS);
			if (p == PosZ.end() ) return(PosZ.size()+1); else return(p-PosZ.begin());
		}

};

const char *DefClusterPath = "C:/Users/Gábor/OneDrive - elte.hu/Documents/Research/pcT/ReadOut/Data/Cluster/";
const char *DefDataPath = "C:/Users/Gábor/OneDrive - elte.hu/Documents/Research/pcT/ReadOut/Data/MC/head_phantom_10000Primaries_230MeV_1617471898.root";

void get_database_root(char const *DataPath = DefDataPath, char const *ClusterPathC = DefClusterPath) {
	std::string ClusterPath = ClusterPathC;
	
	struct Hit {
        float posX;
        float posY;
		float posZ;
        float edep;
		int cS; // clustersize
	};
	std::vector<Hit> Hits;

    ConvertToPixels cP(ClusterPath);

	// read: Data\\MC\\head_phantom_10000Primaries_230MeV_1617471898.root
	// edep, posX, posY
	TFile        * CDB_fCluster;
	cout << "Datafile to process is : " << DataPath << endl;
	CDB_fCluster = new TFile(DataPath, "READ");
	TTreeReader reader("Hits", CDB_fCluster);
	TTreeReaderValue<float> edep(reader, "edep");
	TTreeReaderValue<float> posX(reader, "posX");
	TTreeReaderValue<float> posY(reader, "posY");
	TTreeReaderValue<float> posZ(reader, "posZ");

	Hit t;
	while(reader.Next()) {
		t.posX = *posX;
		t.posY = *posY;
		t.posZ = *posZ;
		t.edep = *edep;
		t.cS = floor(4.2267f*pow(t.edep*40.f,0.65)+0.5f);
		Hits.push_back(t);
	}
	
    //float edep = 0.12f;
//    std::vector<std::array<float,4>> res = cP.get_Pos(0.40f,1.1f,edep);
//    std::vector<std::array<float,4>> res = cP.get_Pos(37.40f,23.50f,edep);
    //std::vector<std::array<float,4>> res = cP.get_Pos(45.02f,-14.60f,edep);
//    std::vector<std::array<float,4>> res = cP.get_Pos(-120.3f,-14.60f,edep);
//    std::vector<std::array<float,4>> res = cP.get_Pos(45.02f,-140.60f,edep);

//	std::size_t ShowMax {10};
//	for(size_t ne=0; ne < min(Hits.size(),ShowMax); ne++) {
	ofstream outfile("nn_samples.txt");
	for(size_t ne=0; ne < Hits.size(); ne++) {
		if(ne < 5) {
			cout << Hits[ne].posZ << " => " << cP.get_Layer(Hits[ne].posZ) << endl;
		}
//		if(cP.get_Layer(Hits[ne].posZ) < 1) {
//			cout << Hits[ne].posZ << " => " << cP.get_Layer(Hits[ne].posZ) << endl;
//		}
		std::vector<std::array<float,4>> res = cP.get_Pos(Hits[ne].posX,Hits[ne].posY,Hits[ne].edep);
		for (size_t i=0; i < res.size(); i++)
		{
			std::copy(res[i].begin(), res[i].end(), 
                std::ostream_iterator<float>(outfile, " "));
			outfile << cP.get_Layer(Hits[ne].posZ) << " " << Hits[ne].edep << " " << Hits[ne].cS << " ";
			std::vector<pair<int,int>> pixs = cP.get_CS_Pos(res[i][2], res[i][3], Hits[ne].edep);
			for (const auto& p : pixs) {
				outfile << p.first << " " << p.second << " ";
			}
			outfile << std::endl;
		}
    }
	outfile.close();

}