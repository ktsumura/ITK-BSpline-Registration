#include <iostream>

#include "BSplineImageRegistrationV4.h"

using namespace std;
using namespace v4;

double * readImage(int fileIndex, unsigned int * imageDim);

int main(int argc, char *argv[])
{
	unsigned int numGridNodes = 10;
	unsigned int numPhases = 20;
	unsigned int imageDim[2] = { 56, 66 };
	unsigned int imageSize = imageDim[0] * imageDim[1];

	double elapsedTime = 0.0;

	for (int phaseIndex = 1; phaseIndex < numPhases; phaseIndex++) {
		// Read the images
		double * fixedImage = readImage(0, imageDim);
		double * movingImage = readImage(phaseIndex, imageDim);

		const clock_t sTime = clock();

		// Create a registration object
		BSplineImageRegistration *registration = new BSplineImageRegistration(fixedImage, movingImage, imageDim);

		// Set the number of grid nodes in on dimension
		registration->SetNumberOfGridNodesInOneDimension(numGridNodes);

		// Register
		registration->Register();

		// Get parameters
		//std::vector<double> parameterVector = registration->GetParameterVector();

		// Dispose
		delete registration;

		const clock_t eTime = clock();
		elapsedTime += (eTime - sTime);
	}
#ifdef WO_VISUALIZATION
	elapsedTime = elapsedTime / CLOCKS_PER_SEC;
	std::cout << "Processing time: " << elapsedTime << std::endl;
#endif 
}

double * readImage(int fileIndex, unsigned int * imageDim) {

	std::string strHeight = std::to_string(imageDim[1]);
	std::string strWidth = std::to_string(imageDim[0]);
	std::string strIndex = std::to_string(fileIndex);
	std::string strFileName = "..\\..\\itkReg\\image\\image_" + strHeight + "_" + strWidth + "_" + strIndex;
	const char * fileName = strFileName.c_str();

	unsigned int size = imageDim[0] * imageDim[1];

	std::ifstream ifs;
	ifs.open(fileName, ios::in | ios::binary);

	double* imageArray = new double[size];

	if (ifs) {
		char * buffer = new char[sizeof(double)];
		char * buffer2 = new char[sizeof(double)];
		for (int i = 0; i < size; i++) {
			// Read 8 bytes
			ifs.read(reinterpret_cast<char*>(buffer), sizeof(double));

			// Swap
			buffer2[0] = buffer[7];
			buffer2[1] = buffer[6];
			buffer2[2] = buffer[5];
			buffer2[3] = buffer[4];
			buffer2[4] = buffer[3];
			buffer2[5] = buffer[2];
			buffer2[6] = buffer[1];
			buffer2[7] = buffer[0];

			imageArray[i] = *reinterpret_cast<double*>(buffer2);
		}
	}
	else {
		std::cout << "Failed to open the file.\n";
	}

	ifs.close();

	return imageArray;
}

