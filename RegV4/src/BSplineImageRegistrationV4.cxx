
#include "BSplineImageRegistrationV4.h"

//#include "itkBSplineTransform.h"

using namespace itk;

namespace v4
{
	BSplineImageRegistration::BSplineImageRegistration(PixelType *fixedArray, PixelType *movingArray, unsigned int* imageDim) {
		m_FixedArrayQueue.push_back(fixedArray);
		m_MovingArrayQueue.push_back(movingArray);

		for (int i = 0; i < ImageDimension; i++) {
			m_ImageDim[i] = imageDim[i];
			m_Origin[i] = 0.0;
			m_PixelSpacing[i] = 1.0;
		}

		m_Metric = ITK_NULLPTR;
		m_Transform = ITK_NULLPTR;
		m_Optimizer = ITK_NULLPTR;
		m_Registration = ITK_NULLPTR;
	}

	BSplineImageRegistration::~BSplineImageRegistration()
	{
		m_FixedArrayQueue.clear();
		m_MovingArrayQueue.clear();
		m_FixedImageQueue.clear();
		m_MovingImageQueue.clear();
		m_Metric = ITK_NULLPTR;
		m_Transform = ITK_NULLPTR;
		m_Optimizer = ITK_NULLPTR;
		m_Registration = ITK_NULLPTR;
	}

	OptimizerType::ParametersType BSplineImageRegistration::GetParameters() {
		return m_Transform->GetParameters();
	}

	CoordinateRepType* BSplineImageRegistration::GetParametersInDouble() {
		TransformType::ParametersType parameters = m_Transform->GetParameters();
		unsigned int numElements = parameters.GetNumberOfElements();

		CoordinateRepType* parameterArray = new CoordinateRepType[numElements];
		for (unsigned int index = 0; index < numElements; index++) {
			parameterArray[index] = parameters.GetElement(index);
		}

		return parameterArray;
	}

	std::vector<CoordinateRepType> BSplineImageRegistration::GetParameterVector() {
		TransformType::ParametersType parameters = m_Transform->GetParameters();
		unsigned int numElements = parameters.GetNumberOfElements();

		std::vector<CoordinateRepType> parameterVector;
		for (unsigned int index = 0; index < numElements; index++) {
			parameterVector.push_back(parameters.GetElement(index));
		}

		return parameterVector;
	}

	void BSplineImageRegistration::AddImages(PixelType *fixedArray, PixelType *movingArray) {
		m_FixedArrayQueue.push_back(fixedArray);
		m_MovingArrayQueue.push_back(movingArray);
	}

	void BSplineImageRegistration::SetNumberOfGridNodesInOneDimension(unsigned int numberOfGridNodesInOneDimension) {
		m_NumberOfGridNodesInOneDimension = numberOfGridNodesInOneDimension;
	}

	void BSplineImageRegistration::Register() {

		// Create fixed images
		if (!m_FixedArrayQueue.empty()) {
			std::deque<PixelType *>::iterator iterator = m_FixedArrayQueue.begin();
			while (iterator != m_FixedArrayQueue.end()) {
				m_FixedImageQueue.push_back(CreateImage(*iterator));
				iterator++;
			}
		}

		// Create moving images
		if (!m_MovingArrayQueue.empty()) {
			std::deque<PixelType *>::iterator iterator = m_MovingArrayQueue.begin();
			while (iterator != m_MovingArrayQueue.end()) {
				m_MovingImageQueue.push_back(CreateImage(*iterator));
				iterator++;
			}
		}

		// Build metric
		BuildMetric();

		// Build scales estimator
		BuildScalesEstimator();

		// Build optimizer
		BuildOptimizer();

		// Build transform
		BuildTransform();

		// Build registration
		BuildRegistration();

		try
		{
			m_Registration->Update();
			//std::cout << "Optimizer stop condition = "
			//	<< m_Registration->GetOptimizer()->GetStopConditionDescription()
			//	<< std::endl;
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return;
		}

		// Get the final parameters
		OptimizerType::ParametersType finalParameters = m_Transform->GetParameters();
#ifdef W_VISUALIZATION
		std::cout << "Last Transform Parameters" << std::endl;
		std::cout << finalParameters << std::endl;

		// The following is for visualization
		ImageType::Pointer fixedImage = m_FixedImageQueue.front();
		ImageType::Pointer movingImage = m_MovingImageQueue.front();

		typedef itk::ResampleImageFilter<ImageType, ImageType > ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetTransform(m_Transform);
		resample->SetInput(movingImage);
		resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
		resample->SetOutputOrigin(fixedImage->GetOrigin());
		resample->SetOutputSpacing(fixedImage->GetSpacing());
		resample->SetOutputDirection(fixedImage->GetDirection());
		resample->SetDefaultPixelValue(100);

		typedef  unsigned char  OutputPixelType;
		typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;

		typedef itk::CastImageFilter<ImageType, OutputImageType > CastFilterType;
		CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput(resample->GetOutput());

		typedef itk::ImageFileWriter< OutputImageType >  OutputWriterType;
		OutputWriterType::Pointer writer = OutputWriterType::New();
		writer->SetFileName("output.png");
		writer->SetInput(caster->GetOutput());

		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return;
		}

		QuickView viewer;
		viewer.AddImage(
			fixedImage.GetPointer(), true,
			"Fixed Image");
		viewer.AddImage(
			movingImage.GetPointer(), true,
			"Moving Image");
		viewer.AddImage(
			resample->GetOutput(), true,
			"Resampled Moving Image");

		viewer.Visualize();
#endif
	}

	ImageType::Pointer BSplineImageRegistration::CreateImage(PixelType *fImage) {
		typedef itk::ImportImageFilter< PixelType, ImageDimension >  ImportFilterType;
		ImportFilterType::Pointer importFilter = ImportFilterType::New();

		// Note imageDim is in the order of row-col
		ImportFilterType::SizeType  size;
		size[0] = m_ImageDim[0];  // size along X
		size[1] = m_ImageDim[1];  // size along Y

		ImportFilterType::IndexType start;
		start.Fill(0);

		ImportFilterType::RegionType region;
		region.SetIndex(start);
		region.SetSize(size);

		importFilter->SetRegion(region);
		importFilter->SetOrigin(m_Origin);
		importFilter->SetSpacing(m_PixelSpacing);

		const itk::SizeValueType imageSize = m_ImageDim[0] * m_ImageDim[1];

		const bool importImageFilterWillOwnTheBuffer = true;
		importFilter->SetImportPointer(fImage, imageSize, importImageFilterWillOwnTheBuffer);
		importFilter->Update();

		ImageType::Pointer  image = ImageType::New();
		image->Graft(importFilter->GetOutput());

		return image;
	}

	void BSplineImageRegistration::BuildMetric() {
		m_Metric = MetricType::New();

		if ((m_FixedImageQueue.size() > 0) && (m_MovingImageQueue.size() > 0)) {
			typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType> MetricCompType;
			MetricCompType::Pointer metricComp1 = MetricCompType::New();
			m_Metric->AddMetric(metricComp1);
		}

		if ((m_FixedImageQueue.size() > 1) && (m_MovingImageQueue.size() > 1)) {
			typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType> MetricCompType;
			MetricCompType::Pointer metricComp2 = MetricCompType::New();
			m_Metric->AddMetric(metricComp2);
		}
	}

	void BSplineImageRegistration::BuildTransform() {
		m_Transform = TransformType::New();
		
		//  Here we define the parameters of the BSplineDeformableTransform grid.  We
		//  arbitrarily decide to use a grid with $5 \times 5$ nodes within the image.
		//  The reader should note that the BSpline computation requires a
		//  finite support region ( 1 grid node at the lower borders and 2
		//  grid nodes at upper borders). Therefore in this example, we set
		//  the grid size to be $8 \times 8$ and place the grid origin such that
		//  grid node (1,1) coincides with the first pixel in the fixed image.
		TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
		for (unsigned int i = 0; i < ImageDimension; i++)
		{
			fixedPhysicalDimensions[i] = m_PixelSpacing[i] * static_cast<CoordinateRepType>(m_ImageDim[i] - 1);
		}
		
		TransformType::MeshSizeType             meshSize;
		meshSize.Fill(m_NumberOfGridNodesInOneDimension - SplineOrder);
		
		m_Transform->SetTransformDomainOrigin(m_Origin);
		m_Transform->SetTransformDomainPhysicalDimensions(fixedPhysicalDimensions);
		m_Transform->SetTransformDomainMeshSize(meshSize);
		
		// Get the number of parameters
		const unsigned int numberOfParameters = m_Transform->GetNumberOfParameters();
		unsigned int numberOfElements = m_Parameters.GetNumberOfElements();
		
		if (numberOfElements != numberOfParameters) {
			// Initialize transform parameters
			//typedef TransformType::ParametersType     ParametersType;
			//ParametersType parameters(numberOfParameters);
			m_Parameters.SetSize(numberOfParameters);
			m_Parameters.Fill(0.0);
		}
				
		// Set transform parameters
		m_Transform->SetParameters(m_Parameters);
	}

	void BSplineImageRegistration::BuildScalesEstimator() {
		m_ScalesEstimator = ScalesEstimatorType::New();
		m_ScalesEstimator->SetMetric(m_Metric);
		m_ScalesEstimator->SetTransformForward(true);
		m_ScalesEstimator->SetSmallParameterVariation(1.0);
	}

	void BSplineImageRegistration::BuildOptimizer() {
		// REGULAR STEP GRADIENT DESCENT
		m_Optimizer = OptimizerType::New();
		//m_Optimizer->SetMaximumStepLength(10.0);
		m_Optimizer->SetMinimumStepLength(0.01);
		m_Optimizer->SetRelaxationFactor(0.7);
		m_Optimizer->SetNumberOfIterations(100);// original 200
		m_Optimizer->SetGradientMagnitudeTolerance(1e-4);
		// Create the Command observer and register it with the optimizer.
		//CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
		//m_Optimizer->AddObserver(itk::IterationEvent(), observer);
	}

	void BSplineImageRegistration::BuildRegistration() {
		m_Registration = RegistrationType::New();

		// The old registration framework has problems with multi-threading
		// For now, we set the number of threads to 1
		m_Registration->SetNumberOfThreads(1);

		m_Registration->SetMetric(m_Metric);
		m_Registration->SetOptimizer(m_Optimizer);
		m_Registration->SetInitialTransform(m_Transform);
		m_Registration->InPlaceOn();

		// Set fixed images
		if (!m_FixedImageQueue.empty()) {
			unsigned int index = 0;
			
			std::deque<ImageType::Pointer>::iterator iterator= m_FixedImageQueue.begin();
			while (iterator != m_FixedImageQueue.end()) {
				m_Registration->SetFixedImage(index, *iterator++);
				index++;
			}
		}

		// Set moving images
		if (!m_MovingImageQueue.empty()) {
			unsigned int index = 0;

			std::deque<ImageType::Pointer>::iterator iterator = m_MovingImageQueue.begin();
			while (iterator != m_MovingImageQueue.end()) {
				m_Registration->SetMovingImage(index, *iterator++);
				index++;
			}
		}

		//  A single level registration process is run using
		//  the shrink factor 1 and smoothing sigma 0.
		//
		const unsigned int numberOfLevels = 1;

		RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
		shrinkFactorsPerLevel.SetSize(1);
		shrinkFactorsPerLevel[0] = 1;

		RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
		smoothingSigmasPerLevel.SetSize(1);
		smoothingSigmasPerLevel[0] = 0;

		m_Registration->SetNumberOfLevels(numberOfLevels);
		m_Registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
		m_Registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	}
}
