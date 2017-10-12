#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <ctime>
#include <algorithm>
#include <ctype.h>
#include <stdexcept>

#define VDW_SF 2.0 //1-4 scaling factor of van der Waals

using namespace std;

static void printWelcomeMessage();
static void showUsage();
static void queryDevices();
bool fileExist(const std::string& name);
__host__ __device__ double boundaryConditions(double sum, const double& dim);
__global__ void calcDerivativesSolventSolvent( double *xArr, double *yArr, double *zArr, const double dimX, const double dimY, const double dimZ, const double watRadius,\
const double watEpsilon, const int waterSite, const int atomSize, const int totAtomSize, const double cut, double *derRadArr, double *derEpsArr, const int numWat); 
static void reportMemStatus();

int main (int argc, char* argv[])
{

    clock_t begin, end;
    double elapsed_secs;
    begin = clock();
    string argvStr, argvStr2;
    string cudaOption = "yes"; 
    string coordsFile;
    string cutoffDistStr;
    double cutoffDist;
    int i, j, k, m, dev;
    int count = 0;
    int numAtoms,  numTotAtoms;
    vector<string> atoms;
    vector<string> atomTypeList;
    vector<double> radList, epsList;
    double watRad, watEps;
    int localBond[2];
    int localAngle[3];
    int localDihed[4];
    string tmp, line, localType;
    double localRad, localEps;

    cudaError_t ierrAsync;
    cudaError_t ierrSync;

    printWelcomeMessage();

    if (argc != 7) {
        showUsage();
        if (argc == 1) {
            cout << "Please provide flags and variables to invoke the program.\n\n"; 
            exit(0); }
        else if (argc == 2){       
            argvStr = argv[1];
            if (argvStr == "-h" || argvStr == "-help" || argvStr == "-Help") 
                exit(0);
            else
                cerr << "Not enough arguments, please try again.\n\n";
                exit(1);}
        else if (argc < 7){ 
            cerr << "Not enough arguments, please try again.\n\n";
            exit(1);}
        else if (argc > 7){
            cerr << "Too many arguments, please try again.\n\n";
            exit(1);}
    } 

    for (i = 1; i < argc; i += 2) 
    {
        if ((i + 1) != argc) {
            argvStr = argv[i]; 
            argvStr2 = argv[i+1];
            transform(argvStr.begin(), argvStr.end(), argvStr.begin(), ::tolower);
            if (argvStr == "-crd")
                coordsFile = argvStr2;
            else if (argvStr == "-cutoff") 
                cutoffDistStr = argvStr2;
            else if (argvStr == "-gpu")
                cudaOption = argvStr2; 
            else {
                showUsage();
                cerr << "Wrong arguments, please try again.\n\n";
                exit(1);}
        }
    }        

    // Check the variables provided in the commandline

    transform(cudaOption.begin(), cudaOption.end(), cudaOption.begin(), ::tolower);

    if (cudaOption.compare("yes")&&cudaOption.compare("no")) {
        cerr << "Aborted. Please use YES or NO for the GPU option." << endl << endl;
        exit(1);
    }  
        
    try {
	cutoffDist = atof(cutoffDistStr.c_str());
    } 
    catch (const std::invalid_argument& e) {
       cerr << "Aborted. Please provide a valid number for the van der Waals cutoff distance." << endl << endl;
       exit(1); }

    if  (cutoffDist < 0) {
        cerr << "Aborted. Please provide a valid number for the van der Waals cutoff distance." << endl << endl;
        exit(1); }

     
    if( !fileExist(coordsFile)) {
        cerr << "Aborted. The trajectory file " << coordsFile << "  does not exist." << endl << endl;
        exit(1);
    } 

    string sumFile = "moleculeInfo";
    string parmFile = "FFparameters.dat";
    
    ifstream sFile(sumFile.c_str());
    ifstream pFile(parmFile.c_str());

    if( !fileExist(sumFile) || !fileExist(parmFile)) {
        cerr << "Aborted. Please use script parseTopoWithParmed.py to generate molecular info first." << endl << endl;
        exit(1);
    }
    else {
        do {
            getline(sFile, line);
        } while (line.find("Number of atoms"));
        
        istringstream ss(line);
        ss >> tmp >> tmp >> tmp >> numTotAtoms;
        cout << "The number of total atoms: " << numTotAtoms << endl;
        
    }
    cout << "Reading in atom types ..." << endl;     

    do {
        getline(sFile, line);
    } while (line.find("[ Atoms ]"));
    
    for (i = 0; i < 4; i++)
        getline(sFile, line);

    do {
        istringstream ss(line);
        ss >> tmp >> tmp >> tmp >> localType;
        atoms.push_back(localType);
        getline(sFile, line);
        count += 1;
    } while (!line.empty());

    numAtoms = count ;
    count = 0;

    cout << "Analyzing connectivity (bonds, angles, 1-4 pairs) of the system ..." << endl;
    cout << "The number of solute atoms (in addition to ions): " << numAtoms << endl;

    // Allocate memory and initialize connectivity matrix
    
    int **connectMatrix = new int*[numAtoms];
    for (i = 0; i < numAtoms; i++) 
        connectMatrix[i] = new int[numAtoms]; 

    for (i = 0; i < numAtoms; i++)
        for (j = 0; j < numAtoms; j++)
            connectMatrix[i][j] = 0;
    
    do {
        getline(sFile, line);
    } while (line.find("[ Bonds ]"));

    getline(sFile, line);
    getline(sFile, line);
    do {
        istringstream ss(line);
        ss >> localBond[0] >> tmp >> tmp >> tmp >> localBond[1];
        connectMatrix[localBond[0]-1][localBond[1]-1] = 1;
        connectMatrix[localBond[1]-1][localBond[0]-1] = 1;
        getline(sFile, line);
    } while (line.find("[ Angles ]"));

    getline(sFile, line);
    getline(sFile, line);
    do {
        istringstream ss(line);
        ss >> localAngle[0] >> tmp >> tmp >> tmp >> localAngle[1] >> tmp >> tmp >> tmp >> localAngle[2];
        connectMatrix[localAngle[0]-1][localAngle[2]-1] = 1;
        connectMatrix[localAngle[2]-1][localAngle[0]-1] = 1;
        getline(sFile, line);
    } while (line.find("[ Dihedrals ]"));
  
    getline(sFile, line);
    getline(sFile, line);
    do {
        istringstream ss(line);
        if ((line[0]!='M')&&(line[0]!='I')) {
            ss >> localDihed[0] >> tmp >> tmp >> tmp >> localDihed[1] >> tmp >> tmp >> tmp >> localDihed[2] >> tmp >> tmp >> tmp >> localDihed[3];
            connectMatrix[localDihed[0]-1][localDihed[3]-1] = 2;
            connectMatrix[localDihed[3]-1][localDihed[0]-1] = 2;
        }
        getline(sFile, line);
    } while (!line.empty());

    cout << "Reading in nonbonded parameters ..." << endl;
    double *radArr = new double[numAtoms];
    double *epsArr = new double[numAtoms];

    do {
        getline(pFile, line);
    } while (line.find("NONB"));
        
    getline(pFile, line);

    do {
        istringstream ss(line);
        ss >> localType >> localRad >> localEps; 
        atomTypeList.push_back(localType);
        radList.push_back(localRad);
        epsList.push_back(localEps);
	getline(pFile, line);
    } while (!line.empty());

    // Assigning nonbonded parameters according to atom types
    for (i = 0; i < numAtoms; i++) {
        for (j = 0; j < atomTypeList.size(); j++) {
            if (!atoms[i].compare(atomTypeList[j])) {
                radArr[i] = radList[j];
                epsArr[i] = epsList[j];
                break;
            }
        }
    }

    for (i = 0; i < atomTypeList.size(); i++) {
        if (!atomTypeList[i].compare("OW")) {
            watRad = radList[i];
            watEps = epsList[i];
            break;
        }
    }
        
    sFile.close();
    pFile.close();    

    /********************  Prepare for the grand loop ***************************/
    cout << "Computing derivatives ..." << endl;
    // reading in coordinates       
    string block;
    int numTotCoordsPerFrame = numTotAtoms * 3;
    int numLines = ceil(numTotCoordsPerFrame/10.0);
    int offset;
    ifstream cFile(coordsFile.c_str());
    cudaError_t ierrDevice = cudaGetDevice( &dev );
    int lineCount;   
    double dim[3];
    int frame = 0;
    double distX, distY, distZ;
    double localDist;
    double radPair, radDist, radDistPow5, radDistPow6, epsPair;
    double rad_der_per_pair, eps_der_per_pair;
    double rad_der_per_type = 0.0;
    double eps_der_per_type = 0.0;
    ofstream rdFile, edFile;
    rdFile.open("radDerivative.dat");
    edFile.open("epsDerivative.dat"); 
    rdFile << " Frame ";
    edFile << " Frame ";
    double wat_rad_per_pair, wat_eps_per_pair;
    double wat_rad_per_type = 0.0;
    double wat_eps_per_type = 0.0;


    // Allocate memory and initialize arrays

    double **coordsMatrix = new double*[numTotAtoms];
    for (i = 0; i < numTotAtoms; i++)
        coordsMatrix[i] = new double[3];

    int atomTypeSize = atomTypeList.size(); 

    string *atomTypeArr = new string[atomTypeSize];
    for (i = 0; i < atomTypeSize; i++){ 
        atomTypeArr[i] = atomTypeList[i];
    }

    for (i = 0; i < atomTypeSize; i++) { 
        rdFile <<  setw(12) << atomTypeList[i];
        edFile <<  setw(12) << atomTypeList[i];
    }
    rdFile << endl;
    edFile << endl;

    string *atomArr = new string[numAtoms];
    for (i = 0; i < numAtoms; i++)
        atomArr[i] = atoms[i];
 
    double *radDerMean = new double[atomTypeSize];
    double *epsDerMean = new double[atomTypeSize];
    for (i = 0; i < atomTypeSize; i++){
        radDerMean[i] = 0;
        epsDerMean[i] = 0;
    }
   
    if (!atomTypeList.empty()) {
        if( atomTypeList.back() == "EP") // Atom type of watervirtual particle found 
            offset = 3;      // for four-site water model
        else 
            offset = 2;    // for three-site water model
    }
         
    
    /********************************* Start the grand loop ********************************************/
    if (!cudaOption.compare("no")) {   
        // CPU code
        getline(cFile, line); 
        getline(cFile, line);
        do {
            frame += 1;
            rdFile << right << setw(7) << frame;
            edFile << right << setw(7) << frame;
            cout << "\r" << "Processing frame " << frame << std::flush; 
            lineCount = 0;
            do {
                lineCount += 1;
                block = block + " " + line;
                getline(cFile, line);
            } while (lineCount!=(numLines + 1));
            istringstream ss(block);
            for (i = 0; i < numTotAtoms; i++) {
                ss >> coordsMatrix[i][0] >> coordsMatrix[i][1] >> coordsMatrix[i][2];
            }    
            ss >> dim[0] >> dim[1] >> dim[2];
          
            for (i = 0; i < atomTypeSize - offset; i++) { // loop over atom types (not including water)
                rad_der_per_type = 0;
                eps_der_per_type = 0;
                
                for (j = 0; j < numAtoms; j++) {
                    if ((!atomTypeArr[i].compare(atomArr[j])) && (epsArr[j] > 0.00001)) {                 
                        for (k = 0; k < numAtoms; k++) { // Compute first-order derivatives involving solute-solute interactions
                            if ((j != k) && (connectMatrix[j][k]!=1) && (epsArr[k] > 0.00001) ) {
                                distX = boundaryConditions(coordsMatrix[j][0] - coordsMatrix[k][0], dim[0]);
                                distY = boundaryConditions(coordsMatrix[j][1] - coordsMatrix[k][1], dim[1]);                                             
                                distZ = boundaryConditions(coordsMatrix[j][2] - coordsMatrix[k][2], dim[2]);
                                localDist = sqrt(distX * distX + distY * distY + distZ * distZ);
                                if (localDist < cutoffDist){
                                    epsPair = sqrt(epsArr[j] * epsArr[k]);
                                    radPair = radArr[j] + radArr[k];
                                    radDist = radPair/localDist;
                                    radDistPow5 = pow(radDist, 5.0);
                                    radDistPow6 = radDistPow5 * radDist;
                                    
                                    rad_der_per_pair = 12.0*(epsPair/localDist)*radDistPow5*(radDistPow6 - 1.0);
                                    eps_der_per_pair = (epsArr[k]/epsPair)*radDistPow6*(0.5*radDistPow6 - 1.0);
 
                                    if (connectMatrix[j][k] == 2) { //1-4 interactions; use a scaling factor of 2 (VDW_SF) for van der Waals                           
                                        rad_der_per_pair = rad_der_per_pair/VDW_SF;
                                        eps_der_per_pair = eps_der_per_pair/VDW_SF;
                                    }
                                    rad_der_per_type += rad_der_per_pair;
                                    eps_der_per_type += eps_der_per_pair;                                                               
                                }
                            }                                     
                        }  // end of the k loop (solute-solute)
    
                        for (k = numAtoms; k < numTotAtoms; k += (offset+1)) { // Considering LJ interactions between solute atoms and water oxygen atoms
                            distX = boundaryConditions(coordsMatrix[j][0] - coordsMatrix[k][0], dim[0]);
                            distY = boundaryConditions(coordsMatrix[j][1] - coordsMatrix[k][1], dim[1]);
                            distZ = boundaryConditions(coordsMatrix[j][2] - coordsMatrix[k][2], dim[2]);
                            localDist = sqrt(distX * distX + distY * distY + distZ * distZ);
                            if (localDist < cutoffDist) {
                                epsPair = sqrt(epsArr[j] * watEps);
                                radPair = radArr[j] + watRad;
                                radDist = radPair/localDist;
                                radDistPow5 = pow(radDist, 5.0);
                                radDistPow6 = radDistPow5 * radDist;

                                rad_der_per_pair = 12.0*(epsPair/localDist)*radDistPow5*(radDistPow6 - 1.0);
                                eps_der_per_pair = (watEps/epsPair)*radDistPow6*(0.5*radDistPow6 - 1.0);
                                
                                wat_rad_per_pair = rad_der_per_pair;
                                wat_eps_per_pair = eps_der_per_pair * epsArr[j] / watEps; 

                                rad_der_per_type += rad_der_per_pair;
                                eps_der_per_type += eps_der_per_pair;
                                wat_rad_per_type += wat_rad_per_pair;
                                wat_eps_per_type += wat_eps_per_pair;

  
                            }
                        } // end of the k loop (solute-solvent)
                    }
                              
                } // end of j loop
                rdFile << right << setw(12) << setprecision(4) << fixed << rad_der_per_type;
                edFile << right << setw(12) << setprecision(4) << fixed << eps_der_per_type;
   
                radDerMean[i] = radDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + rad_der_per_type / (frame * 1.0);
                epsDerMean[i] = epsDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + eps_der_per_type / (frame * 1.0);

            } // end of i loop

            for (j = numAtoms; j < numTotAtoms; j += (offset+1)) {
                for (k = j + offset + 1; k < numTotAtoms; k += (offset+1)) {
                    distX = boundaryConditions(coordsMatrix[j][0] - coordsMatrix[k][0], dim[0]);
                    distY = boundaryConditions(coordsMatrix[j][1] - coordsMatrix[k][1], dim[1]);
                    distZ = boundaryConditions(coordsMatrix[j][2] - coordsMatrix[k][2], dim[2]);
                    localDist = sqrt(distX * distX + distY * distY + distZ * distZ);
                    if (localDist < cutoffDist){
                        radDist = watRad * 2.0/localDist;
                        radDistPow5 = pow(radDist, 5.0);
                        radDistPow6 = radDistPow5 * radDist;
                        wat_rad_per_pair = 12.0*(watEps/localDist)*radDistPow5*(radDistPow6 - 1.0);
                        wat_eps_per_pair = radDistPow6*(0.5*radDistPow6 - 1.0);

                        wat_rad_per_type += wat_rad_per_pair * 2.0;
                        wat_eps_per_type += wat_eps_per_pair * 2.0;
                    } 
                }
            }

            rdFile << right << setw(12) << setprecision(4) << fixed << wat_rad_per_type;
            edFile << right << setw(12) << setprecision(4) << fixed << wat_eps_per_type;
                       
            radDerMean[i] = radDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + wat_rad_per_type / (frame * 1.0);
            epsDerMean[i] = epsDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + wat_eps_per_type / (frame * 1.0);

            wat_rad_per_type = 0.0;
            wat_eps_per_type = 0.0;

            block.clear();
            rdFile << endl;
            edFile << endl;
        } while (!cFile.eof());
        
        rdFile << "   Mean";
        edFile << "   Mean";
        
        for (m = 0; m < atomTypeSize; m++) {
            rdFile <<  right << setw(12) << setprecision(4) << fixed << radDerMean[m];
            edFile <<  right << setw(12) << setprecision(4) << fixed << epsDerMean[m];
        }

        rdFile << endl;
        edFile << endl;

        cout << endl;
    }
    else {
        // CUDA CODE
        queryDevices();    
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, dev);
        
        double *xCoords, *yCoords, *zCoords;
        double *derArrEps, *derArrRad;    

        // Allocate memory
        int cudaArrSize = numTotAtoms *  sizeof (double);
        cudaMallocManaged(&xCoords, cudaArrSize);
        cudaMallocManaged(&yCoords, cudaArrSize);
        cudaMallocManaged(&zCoords, cudaArrSize);

        int nx = (numTotAtoms - numAtoms)/(offset + 1);
        int ny = nx;

        int derArrSize = nx * ny * sizeof (double);
        cudaMallocManaged(&derArrRad, derArrSize);
        cudaMallocManaged(&derArrEps, derArrSize);

        reportMemStatus();

        for (m = 0; m < nx * ny; m++){
            derArrRad[m] = 0;
            derArrEps[m] = 0;
        }

        dim3 threads_per_block (16, 16, 1);
        dim3 number_of_blocks ((nx + threads_per_block.x -1)/ threads_per_block.x, (ny + threads_per_block.y - 1)/ threads_per_block.y, 1);

        getline(cFile, line);
        getline(cFile, line);
        do {
            frame += 1;
            rdFile << right << setw(7) << frame;
            edFile << right << setw(7) << frame;
            cout << "\r" << "Processing frame " << frame << std::flush;
            lineCount = 0;
            do {
                lineCount += 1;
                block = block + " " + line;
                getline(cFile, line);
            } while (lineCount!=(numLines + 1));
            istringstream ss(block);
            for (i = 0; i < numTotAtoms; i++) {
               	ss >> xCoords[i] >> yCoords[i] >> zCoords[i];
            }

            ss >> dim[0] >> dim[1] >> dim[2];
            for (i = 0; i < atomTypeSize - offset; i++) { // loop over atom types (not including water)
                rad_der_per_type = 0;
                eps_der_per_type = 0;

                for (j = 0; j < numAtoms; j++) {
                    if ((!atomTypeArr[i].compare(atomArr[j])) && (epsArr[j] > 0.00001)) {
                        for (k = 0; k < numAtoms; k++) { // Compute first-order derivatives involving solute-solute interactions
                            if ((j != k) && (connectMatrix[j][k]!=1) && (epsArr[k] > 0.00001) ) {
                                distX = boundaryConditions(xCoords[j] - xCoords[k], dim[0]);
                                distY = boundaryConditions(yCoords[j] - yCoords[k], dim[1]);
                                distZ = boundaryConditions(zCoords[j] - zCoords[k], dim[2]);
                                localDist = sqrt(distX * distX + distY * distY + distZ * distZ);
                                if (localDist < cutoffDist){
                                    epsPair = sqrt(epsArr[j] * epsArr[k]);
                                    radPair = radArr[j] + radArr[k];
                                    radDist = radPair/localDist;
                                    radDistPow5 = pow(radDist, 5.0);
                                    radDistPow6 = radDistPow5 * radDist;
                                    rad_der_per_pair = 12.0*(epsPair/localDist)*radDistPow5*(radDistPow6 - 1.0);
                                    eps_der_per_pair = (epsArr[k]/epsPair)*radDistPow6*(0.5*radDistPow6 - 1.0);

                                    if (connectMatrix[j][k] == 2) { //1-4 interactions; use a scaling factor of 2 (VDW_SF) for van der Waals
                                        rad_der_per_pair = rad_der_per_pair/VDW_SF;
                                        eps_der_per_pair = eps_der_per_pair/VDW_SF;
                                    }
                                    rad_der_per_type += rad_der_per_pair;
                                    eps_der_per_type += eps_der_per_pair;
                                }
                            }
                        }  // end of the k loop (solute-solute)

                        for (k = numAtoms; k < numTotAtoms; k += (offset+1)) { // Considering LJ interactions between solute atoms and water oxygen atoms
                            distX = boundaryConditions(xCoords[j] - xCoords[k], dim[0]);
                            distY = boundaryConditions(yCoords[j] - yCoords[k], dim[1]);
                            distZ = boundaryConditions(zCoords[j] - zCoords[k], dim[2]);
                            localDist = sqrt(distX * distX + distY * distY + distZ * distZ);
                            if (localDist < cutoffDist) {
                                epsPair = sqrt(epsArr[j] * watEps);
                                radPair = radArr[j] + watRad;
                                radDist = radPair/localDist;
                                radDistPow5 = pow(radDist, 5.0);
                                radDistPow6 = radDistPow5 * radDist;

                                rad_der_per_pair = 12.0*(epsPair/localDist)*radDistPow5*(radDistPow6 - 1.0);
                                eps_der_per_pair = (watEps/epsPair)*radDistPow6*(0.5*radDistPow6 - 1.0);

                                wat_rad_per_pair = rad_der_per_pair;
                                wat_eps_per_pair = eps_der_per_pair * epsArr[j] / watEps;                                
 
                                rad_der_per_type += rad_der_per_pair;
                                eps_der_per_type += eps_der_per_pair;
                                wat_rad_per_type += wat_rad_per_pair;
                                wat_eps_per_type += wat_eps_per_pair;

                             }
                        } // end of the k loop (solute-solvent)
                    }

                } // end of j loop
                rdFile << right << setw(12) << setprecision(4) << fixed << rad_der_per_type;
                edFile << right << setw(12) << setprecision(4) << fixed << eps_der_per_type;

                radDerMean[i] = radDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + rad_der_per_type / (frame * 1.0);
                epsDerMean[i] = epsDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + eps_der_per_type / (frame * 1.0);
            } // end of i loop

            calcDerivativesSolventSolvent <<< number_of_blocks, threads_per_block >>> (xCoords, yCoords, zCoords, dim[0], dim[1], dim[2], watRad, watEps, \
                                     offset+1, numAtoms, numTotAtoms, cutoffDist, derArrRad, derArrEps, nx);

            ierrSync = cudaGetLastError();
            if ( cudaSuccess != ierrSync ) {
                printf( "Kernel did not launch!\n" );
                printf("Sync error: %s\n", cudaGetErrorString(ierrSync)); 
                exit(1);
             }        

            ierrAsync = cudaDeviceSynchronize(); // Wait for the GPU to finish

            if (ierrAsync != cudaSuccess) { printf("Async error: %s\n", cudaGetErrorString(ierrAsync)); exit(1);} 

            for (m = 0; m < nx*ny; m++) {

                wat_rad_per_type += derArrRad[m];
                wat_eps_per_type += derArrEps[m];

            }
            rdFile << right << setw(12) << setprecision(4) << fixed << wat_rad_per_type;
            edFile << right << setw(12) << setprecision(4) << fixed << wat_eps_per_type;
                       
            radDerMean[i] = radDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + wat_rad_per_type / (frame * 1.0);
            epsDerMean[i] = epsDerMean[i] * ((frame - 1)*1.0) / (frame * 1.0) + wat_eps_per_type / (frame * 1.0);

            wat_rad_per_type = 0.0;
            wat_eps_per_type = 0.0;
             
            block.clear();
    
            rdFile << endl;
            edFile << endl;
        } while (!cFile.eof());
        
        // Free all our allocated memory
        cudaFree( xCoords ); cudaFree( yCoords ); cudaFree( zCoords );
        cudaFree(derArrRad);
        cudaFree(derArrEps);

        rdFile << "   Mean";
        edFile << "   Mean";

        for (m = 0; m < atomTypeSize; m++) {
            rdFile <<  right << setw(12) << setprecision(4) << fixed << radDerMean[m];
            edFile <<  right << setw(12) << setprecision(4) << fixed << epsDerMean[m];
        }

        rdFile << endl;
        edFile << endl;

        cout << endl;
    }
    // Clean up 
    for (i = 0; i < numAtoms; i++)
        delete [] connectMatrix[i];

    for (i = 0; i < numTotAtoms; i++)
        delete [] coordsMatrix[i];
    
    delete [] connectMatrix;
    delete [] coordsMatrix;
    delete [] radArr;
    delete [] epsArr;
    delete [] atomTypeArr;
    delete [] atomArr;
    delete [] radDerMean;
    delete [] epsDerMean;

    rdFile.close();
    edFile.close();

    end = clock();
    elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
    double elapsed_hours = elapsed_secs / 3600.0;
    cout << "It took " << setprecision(2) << fixed << elapsed_hours << " hours to compute the nonbonded LJ derivatives." << endl;

}    

static void printWelcomeMessage()
{
    std::cout << "**********************************************************************\n"
              << "                 Sensitivity Analysis            \n"
              << "Written by:        Jian (Jane) Yin               \n"
              << "under the direction of Prof. Michael K. Gilson   \n"
              << "                     September 2017              \n"
              << "**********************************************************************\n"
              << endl;
}


bool fileExist (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}


static void showUsage()
{
    std::cout << "Flags:\n"
              << "  -crd       MD trajectory file\n"
              << "  -cutoff    cutoff distance of nonbonded interactions\n"
              << "  -gpu       yes or no\n"
              <<"\n"
              << "Example:\n"
              << "  ./sensitivity -crd traj.mdcrd -cutoff 9.0 -gpu yes\n"
              << endl;
}

__host__ __device__ double boundaryConditions(double sum, const double& dim)
{
    if (sum <= -dim*0.5)
        sum += dim;
    else if (sum > dim*0.5)
        sum -= dim;
    return sum;
}

static void queryDevices(){
    // Query GPU information
    int driverVersion = 0;
    int runtimeVersion = 0;
    int dev;
    int deviceCount = 0;
    cudaDeviceProp prop;

    cout << "Sensitivity_CUDA was enabled." << endl;

    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess) {
        printf("%s.\n",cudaGetErrorString(error_id));
        cout << "Aborted." << endl << endl;
        exit(1);
    }

    if (deviceCount == 0) {
        cout << endl << "Aborted. There are no available device(s) that support CUDA." << endl << endl;
        exit(1);
    }

    cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);
    printf("The CUDA Driver Version / Runtime Version is %d / %d.\n", (int)(driverVersion/1000.0),  (int)(runtimeVersion/1000.0));

    cudaGetDevice( &dev );
    cudaGetDeviceProperties(&prop, dev);
    printf("Device %d: \"%s\".\n", dev, prop.name);
    printf("Compute capability: %d.%d.\n", prop.major, prop.minor);

    cout << endl;
}


__global__ void calcDerivativesSolventSolvent( double *xArr, double *yArr, double *zArr, const double dimX, const double dimY, const double dimZ, const double watRadius,\
const double watEpsilon, const int waterSite, const int atomSize, const int totAtomSize, const double cut, double *derRadArr, double *derEpsArr, const int numWat) 
{              
           
   
    double distX, distY, distZ;
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    double localDist, radDist, radDistPow5, radDistPow6;

    int xc, yc;

    if (gx < numWat && gy < numWat) {
        
        if (gx != gy) {

            xc = atomSize + gx * waterSite;
            yc = atomSize + gy * waterSite; 
        
            distX = boundaryConditions(xArr[xc] - xArr[yc], dimX);
            distY = boundaryConditions(yArr[xc] - yArr[yc], dimY);
            distZ = boundaryConditions(zArr[xc] - zArr[yc], dimZ);
            localDist = sqrt(distX * distX + distY * distY + distZ * distZ);
            
            if (localDist < cut){
                radDist = watRadius * 2.0/localDist;
                radDistPow5 = pow(radDist, 5.0);
                radDistPow6 = radDistPow5 * radDist;
               
                derRadArr[gx*numWat + gy] = 12.0*(watEpsilon/localDist)*radDistPow5*(radDistPow6 - 1.0);
                derEpsArr[gx*numWat + gy] = radDistPow6*(0.5*radDistPow6 - 1.0);
            }
            else {
                derRadArr[gx*numWat + gy] = 0.0;
               	derEpsArr[gx*numWat + gy] = 0.0;
           }  
                

        }
        
    }
   
}

/**
 * API to report the memory usage of the GPU
 */
static void reportMemStatus() {

	// show memory usage of GPU
	size_t free_byte;
	size_t total_byte;
	size_t malloc_byte;

	cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

	if (cudaSuccess != cuda_status) {
		printf("Error: cudaMemGetInfo fails, %s \n",
				cudaGetErrorString(cuda_status));
		return;
	}

	cuda_status = cudaDeviceGetLimit(&malloc_byte, cudaLimitMallocHeapSize);
	if (cudaSuccess != cuda_status) {
			printf("Error: cudaDeviceGetLimit fails, %s \n",
					cudaGetErrorString(cuda_status));
			return;
	}

	double free_db = (double) free_byte;
	double total_db = (double) total_byte;
	double used_db = total_db - free_db;
	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB, malloc limit = %f MB\n",
			used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0,
			total_db / 1024.0 / 1024.0, malloc_byte / 1024.0 / 1024.0);

}

