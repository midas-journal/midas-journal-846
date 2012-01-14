#include "itkPhaseSymmetryImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkArray2D.h"
#include <vector>
#include <string>
#include <sstream>
#include <ctime>
#include <exception>

std::vector<std::string> parse(std::string l, char delim) 
{
	std::replace(l.begin(), l.end(), delim, ' ');
	std::istringstream stm(l);
	std::vector<std::string> tokens;
	for (;;) {
		std::string word;
		if (!(stm >> word)) break;
		tokens.push_back(word);
	}
	return tokens;
}

int main(int argc, char* argv[])
{

	std::stringstream ss;
	std::string infile;
	std::string outfile;
	std::string wavelengthsS;
	std::string orientationsS;
	std::string sigmaS;
	std::string angular_bandwidthS;
	std::string polarityS;
	std::string ntS;
	
	double pi = 3.1416;
/*
	if (argc < 3)
	{
		std::cerr << "Usage: PhaseSymmetryFilter3D.exe infile outfile wavelengths orientations sigma angular_bandwidth polarity noise_threshhold" << std::endl;
		std::cerr << "Example: PhaseSymmetryFilter3D.exe i.mhd o.mhd 3,3,3,6,6,6,12,12,12 1,0,0,0,1,0,0,0,1 0.55 3.14 0 10.0" << std::endl;
		return EXIT_FAILURE;
	}*/
	if(argc < 3)
	{
		ss << "C:\\data\\test3D.mhd";
		infile = ss.str();
		ss.str("");
		ss << "C:\\data\\test3DOut.mhd";
		outfile = ss.str();
		ss.str("");
		ss << "3,3,3,6,6,6,12,12,12";
		wavelengthsS = ss.str();
		ss.str("");
		ss << "1,0,0,0,1,0,0,0,1";
		orientationsS = ss.str();
		ss.str("");
		ss << "0.55";
		sigmaS = ss.str();
		ss.str("");
		ss << "3.1416";
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << "1";
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else if(argc < 4)
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << "3,3,3,6,6,6,12,12,12";
		wavelengthsS = ss.str();
		ss.str("");
		ss << "1,0,0,0,1,0,0,0,1";
		orientationsS = ss.str();
		ss.str("");
		ss << "0.55";
		sigmaS = ss.str();
		ss.str("");
		ss << "3.1416";
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << "1";
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else if(argc < 5)
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << argv[3];
		wavelengthsS = ss.str();
		ss.str("");
		ss << "1,0,0,0,1,0,0,0,1";
		orientationsS = ss.str();
		ss.str("");
		ss << "0.55";
		sigmaS = ss.str();
		ss.str("");
		ss << "3.1416";
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << "1";
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else if(argc < 6)
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << argv[3];
		wavelengthsS = ss.str();
		ss.str("");
		ss << argv[4];
		orientationsS = ss.str();
		ss.str("");
		ss << "0.55";
		sigmaS = ss.str();
		ss.str("");
		ss << "3.1416";
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << "1";
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else if(argc < 7)
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << argv[3];
		wavelengthsS = ss.str();
		ss.str("");
		ss << argv[4];
		orientationsS = ss.str();
		ss.str("");
		ss << argv[5];
		sigmaS = ss.str();
		ss.str("");
		ss << "3.1416";
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << "1";
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else if(argc < 8)
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << argv[3];
		wavelengthsS = ss.str();
		ss.str("");
		ss << argv[4];
		orientationsS = ss.str();
		ss.str("");
		ss << argv[5];
		sigmaS = ss.str();
		ss.str("");
		ss << argv[6];
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << "1";
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else if(argc < 9)
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << argv[3];
		wavelengthsS = ss.str();
		ss.str("");
		ss << argv[4];
		orientationsS = ss.str();
		ss.str("");
		ss << argv[5];
		sigmaS = ss.str();
		ss.str("");
		ss << argv[6];
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << argv[7];
		polarityS = ss.str();
		ss.str("");
		ss << "10.0";
		ntS = ss.str();
		ss.str("");
	}
	else
	{
		ss << argv[1];
		infile = ss.str();
		ss.str("");
		ss << argv[2];
		outfile = ss.str();
		ss.str("");
		ss << argv[3];
		wavelengthsS = ss.str();
		ss.str("");
		ss << argv[4];
		orientationsS = ss.str();
		ss.str("");
		ss << argv[5];
		sigmaS = ss.str();
		ss.str("");
		ss << argv[6];
		angular_bandwidthS = ss.str();
		ss.str("");
		ss << argv[7];
		polarityS = ss.str();
		ss.str("");
		ss << argv[8];
		ntS = ss.str();
		ss.str("");
	}


	typedef float ImagePixelType;
	int ndims = 3;
	typedef itk::Image< ImagePixelType, 3 > ImageType;
	typedef itk::PhaseSymmetryImageFilter< ImageType, ImageType > PSFilterType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;

	typedef itk::Array2D<double> DoubleMatrix;

	ImageType::Pointer inImg = ImageType::New();
	ReaderType::Pointer reader  = ReaderType::New();
	WriterType::Pointer writer= WriterType::New();
	PSFilterType::Pointer psfilter = PSFilterType::New();


	//3 wavelengths, 3 dimensions
	std::vector<std::string> wavelengthsVS = parse(wavelengthsS,',');
	if(wavelengthsVS.size()%3!=0)
	{
		std::cerr << "wavelengths must be a comma seperated string of numbers with number of elements divisible by 3" << std::endl;
		return EXIT_FAILURE;
	}

	std::vector<std::string> orientationsVS = parse(orientationsS,',');
	if(orientationsVS.size()%3!=0)
	{
		std::cerr << "orientations must be a comma seperated string of numbers with number of elements divisible by 3" << std::endl;
		return EXIT_FAILURE;
	}

	int wvCount = int(double(wavelengthsVS.size())/3.0); 
	DoubleMatrix wavelengths(wvCount,3);
	int orCount = int(double(orientationsVS.size())/3.0); 
	DoubleMatrix orientations(orCount,3);

	int idx=0;
	for(int i=0; i < wvCount; i++)
	{
		for(int j=0; j < 3; j++)
		{
			wavelengths(i,j) = atof(wavelengthsVS[idx].c_str());
			idx++;
		}
	}
	

	idx=0;
	for(int i=0; i < orCount; i++)
	{
		for(int j=0; j < 3; j++)
		{
			orientations(i,j) = atof(orientationsVS[idx].c_str());
			idx++;
		}
	}
	

	double sigma=0;
	sigma = atof( sigmaS.c_str() ) ;


	double anglebandwidth=0;
	anglebandwidth = atof( angular_bandwidthS.c_str() ) ;


	int polarity=0;
	polarity =int( atof( polarityS.c_str() ) );


	double noiseT=0;
	noiseT = atof( ntS.c_str() ) ;

	//std::cerr << infile.c_str() << std::endl;
	//std::cerr << outfile.c_str() << std::endl;

	reader->SetFileName(infile.c_str());
	try
	{
		reader->Update();
		inImg = reader->GetOutput();
		inImg->DisconnectPipeline();
	}
	catch ( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	psfilter->SetInput(inImg);
	psfilter->SetWavelengths(wavelengths);
	psfilter->SetOrientations(orientations);
	psfilter->SetSigma(sigma);
	psfilter->SetAngleBandwidth(anglebandwidth);
	psfilter->SetPolarity(polarity);
	psfilter->SetT(noiseT);
	psfilter->Initialize();


	try
	{
		writer->SetFileName(outfile.c_str());
		writer->SetInput(psfilter->GetOutput());
		writer->Update();
	}
	catch ( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}


