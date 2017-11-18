#ifndef _BSPLINE_IMAGE_REGISTRATION_
#define _BSPLINE_IMAGE_REGISTRATION_

#include <deque>;

#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkImportImageFilter.h"
#include "QuickView.h"
#include "itkTransformToDisplacementFieldFilter.h"

using namespace std;

namespace v3
{
	const unsigned int ImageDimension = 2;
	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;

	typedef double PixelType;
	typedef double CoordinateRepType;

	typedef itk::Image< PixelType, ImageDimension >  ImageType;
	typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
	typedef itk::BSplineTransform<CoordinateRepType, SpaceDimension, SplineOrder > TransformType;
	typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateRepType>    InterpolatorType;
	typedef itk::ImageRegistrationMethod<ImageType, ImageType >    RegistrationType;

	class BSplineImageRegistration
	{
	public:
		BSplineImageRegistration(PixelType *m_FixedArray, PixelType *m_MovingArray, unsigned int* m_imageDim);
		~BSplineImageRegistration();

		OptimizerType::ParametersType GetParameters();
		CoordinateRepType* GetParametersInDouble();
		std::vector<CoordinateRepType> GetParameterVector();

		void AddImages(PixelType *fixedArray, PixelType *movingArray);
		void SetNumberOfGridNodesInOneDimension(unsigned int numberOfGridNodesInOneDimension);

		void Register();
	protected:
	private:
		std::deque<PixelType *>            m_FixedArrayQueue;
		std::deque<PixelType *>            m_MovingArrayQueue;

		itk::SpacePrecisionType m_ImageDim[ImageDimension];
		itk::SpacePrecisionType m_Origin[ImageDimension];
		itk::SpacePrecisionType m_PixelSpacing[ImageDimension];
		TransformType::ParametersType m_Parameters;

		std::map<std::string, std::vector<CoordinateRepType>> m_FixedContourX;
		std::map<std::string, std::vector<CoordinateRepType>> m_FixedContourY;
		std::map<std::string, std::vector<CoordinateRepType>> m_MovingContourX;
		std::map<std::string, std::vector<CoordinateRepType>> m_MovingContourY;

		unsigned int m_NumberOfGridNodesInOneDimension;

		std::deque<ImageType::Pointer>            m_FixedImageQueue;
		std::deque<ImageType::Pointer>            m_MovingImageQueue;

		MetricType::Pointer m_Metric;
		TransformType::Pointer  m_Transform;
		OptimizerType::Pointer m_Optimizer;
		InterpolatorType::Pointer m_Interpolator;
		RegistrationType::Pointer m_Registration;

		void BuildMetric();
		void BuildTransform();
		void BuildOptimizer();
		void BuildInterpolator();
		void BuildRegistration();
		ImageType::Pointer CreateImage(PixelType *fImage);


		class CommandIterationUpdate : public itk::Command
		{
		public:
			typedef  CommandIterationUpdate   Self;
			typedef  itk::Command             Superclass;
			typedef itk::SmartPointer<Self>   Pointer;
			itkNewMacro(Self);

		protected:
			CommandIterationUpdate() {};

		public:
			typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
			typedef   const OptimizerType *                  OptimizerPointer;

			void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
			{
				Execute((const itk::Object *)caller, event);
			}

			void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
			{
				OptimizerPointer optimizer =
					static_cast< OptimizerPointer >(object);
				if (!(itk::IterationEvent().CheckEvent(&event)))
				{
					return;
				}
				std::cout << "Iteration : ";
				std::cout << optimizer->GetCurrentIteration() << "   ";
				std::cout << optimizer->GetValue() << "   ";
				std::cout << std::endl;
			}
		};
	};
}

#endif _BSPLINE_IMAGE_REGISTRATION_