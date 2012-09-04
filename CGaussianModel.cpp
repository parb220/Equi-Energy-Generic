#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <cstring>
#include "CGaussianModel.h"

void EigenAnalysis(double *, double *, int, const double *);

CGaussianModel::CGaussianModel(int _dim, const double *_mean, const double *_covariance) : CSimpleGaussianModel(_dim, _mean, NULL)
{
	if (nData > 0)
	{
		covariance_matrix = new double[nData*nData]; 
		transformation_matrix = new double[nData*nData]; 
		if (_covariance != NULL)
		{
			memcpy(covariance_matrix, _covariance, nData*nData*sizeof(double)); 
			EigenAnalysis(transformation_matrix, sigma, nData, covariance_matrix);
			TransformMean(); 
		}
		else // by default, covariance_matrix = I
		{
			for (int i=0; i<nData; i++)
			{
				for (int j=0; j<nData && j!=i; j++)
					covariance_matrix[i*nData+j] = transformation_matrix[i*nData+j] = 0.0;
				covariance_matrix[i*nData+i] = transformation_matrix[i*nData+i] = 1.0;  
			}
		}
	}
}

CGaussianModel::CGaussianModel(const CGaussianModel &right) : CSimpleGaussianModel(right.nData, right.mu, right.sigma)
{
	if (nData > 0)
	{
		covariance_matrix = new double[nData*nData]; 
		transformation_matrix = new double[nData*nData]; 
		memcpy(covariance_matrix, right.covariance_matrix, nData*nData*sizeof(double)); 
		memcpy(transformation_matrix, right.transformation_matrix, nData*nData*sizeof(double)); 
	}
}

const CGaussianModel & CGaussianModel::operator = (const CGaussianModel &right)
{
	SetDataDimension(right.nData); 
	SetMeanParameter(right.mu, nData); 
	SetSigmaParameter(right.sigma, nData); 
	memcpy(covariance_matrix, right.covariance_matrix, nData*nData*sizeof(double));
        memcpy(transformation_matrix, right.transformation_matrix, nData*nData*sizeof(double));
	return *this; 
}

CGaussianModel::~CGaussianModel()
{
	if (nData >0)
	{
		delete [] covariance_matrix; 
		delete [] transformation_matrix; 
	}
}

void CGaussianModel::SetDataDimension(int _dim)
{
	if (nData != _dim)
	{
		if (nData > 0)
		{
			delete [] covariance_matrix; 
			delete [] transformation_matrix;
		}
		CSimpleGaussianModel::SetDataDimension(_dim); 
		covariance_matrix = new double[nData*nData];
                transformation_matrix = new double[nData*nData];
		for (int i=0; i<nData; i++)
                {
			for (int j=0; j<nData && j!=i; j++)
				covariance_matrix[i*nData+j] = transformation_matrix[i*nData+j] = 0.0;
			covariance_matrix[i*nData+i] = transformation_matrix[i*nData+i] = 1.0;
                }
	}
}

void CGaussianModel::SetCovarianceMatrix(const double *_covariance, int _dim)
{
	SetDataDimension(_dim); 
	memcpy(covariance_matrix, _covariance, nData*nData*sizeof(double)); 
	EigenAnalysis(transformation_matrix, sigma, nData, covariance_matrix);
	TransformMean(); 
}

void CGaussianModel::GetCovarianceMatrix(double *_buffer, int _dim) const
{
	memcpy(_buffer, covariance_matrix, nData*nData*sizeof(double)); 
}

void CGaussianModel::SetTransformationMatrix(const double *_transformation, int _dim)
{
	SetDataDimension(_dim); 
	memcpy(transformation_matrix, _transformation, nData*nData*sizeof(double)); 
	TransformMean(); 
}

void CGaussianModel::GetEigenVector(double *_eigenV, int nDim, int jIndex) const 
{
	// transformation_matrix(:, jIndex)
	for (int i=0; i<nData; i++)
		_eigenV[i] = transformation_matrix[i*nData+jIndex];  
}

void CGaussianModel::GetEigenMatrix_Row(double *_row, int nDim, int iIndex) const 
{
	// transformation_matrix(iIndex, :)
	memcpy(_row, transformation_matrix+iIndex*nData, nData*sizeof(double)); 
}

void EigenAnalysis(double *eVector, double *eValue, int dim, const double *covariance)
{
	double *covariance_view = new double[dim*dim]; 
	memcpy(covariance_view, covariance, dim*dim*sizeof(double)); 
	gsl_matrix_view m = gsl_matrix_view_array(covariance_view, dim, dim); 
	gsl_vector *eval = gsl_vector_alloc(dim); 
	gsl_matrix *evec = gsl_matrix_alloc(dim, dim); 
	gsl_eigen_symmv_workspace *w= gsl_eigen_symmv_alloc(dim); 
	gsl_eigen_symmv(&m.matrix, eval, evec, w); 
	gsl_eigen_symmv_free(w); 
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC); 
	for (int i=0; i<dim; i++)
	{
		eValue[i] = gsl_vector_get(eval, i); 
		for (int j=0; j<dim; j++)
			eVector[i*dim+j] = gsl_matrix_get(evec, i, j); 
	}
	delete []covariance_view; 
	gsl_vector_free(eval); 
	gsl_matrix_free(evec); 
}

double InnerProduct(const double *x, const double *y, int dim)
{
	double result = 0; 
	for (int i=0; i<dim; i++)
		result += x[i]*y[i]; 
	return result; 
}

void CGaussianModel::RotateByTransformationMatrix(double *y, int dY, const double *x) const
{
	double *eVector = new double [nData]; 
	for (int i=0; i<nData; i++)
	{
		GetEigenVector(eVector, nData, i); 
		y[i] = InnerProduct(eVector, x, nData); 
	}
	delete [] eVector; 
}

void CGaussianModel::InverseRotateByTransformationMatrix(double *y, int dY, const double *x) const
{
	double *rowVector = new double [nData]; 
	for (int i=0; i<nData; i++)
	{
		GetEigenMatrix_Row(rowVector, nData, i); 
		y[i] = InnerProduct(rowVector, x, nData); 
	}
	delete [] rowVector; 
}

void CGaussianModel::TransformMean()
{
	double *mu_rotated = new double [nData]; 
	RotateByTransformationMatrix(mu_rotated, nData, mu); 
	SetMeanParameter(mu_rotated, nData); 
	delete [] mu_rotated; 
}

double CGaussianModel::log_prob(CSampleIDWeight &x) const
{
	CSampleIDWeight x_rotated; 
	x_rotated.SetDataDimension(nData);   
	RotateByTransformationMatrix(x_rotated.GetData(), nData, x.GetData()); 
	double result = CSimpleGaussianModel::log_prob(x_rotated); 
	x.log_prob = x_rotated.log_prob; 
	x.SetWeight(x_rotated.GetWeight()); 
	return result; 
}

void CGaussianModel::draw(CSampleIDWeight &y, bool &new_sample_flag, const gsl_rng* r, int mMH) const
{
	CSampleIDWeight y_rotated; 
	y_rotated.SetDataDimension(nData); 
	CSimpleGaussianModel::draw(y_rotated, new_sample_flag, r); 

	y.SetDataDimension(nData); 
	InverseRotateByTransformationMatrix(y.GetData(), nData, y_rotated.GetData()); 
	y.log_prob = y_rotated.log_prob; 
	y.SetWeight(y_rotated.GetWeight()); 
}

void CGaussianModel::GetMode(CSampleIDWeight &x, int iMode) const
{
	CSampleIDWeight x_rotated(mu,nData,0, 0.0);  
	// mu has been rotated before and has to be rotated back now
	CSimpleGaussianModel::log_prob(x_rotated);	// directly calls SimpleGaussian::log_prob to calculate log_prob  
	InverseRotateByTransformationMatrix(x.GetData(), nData, x_rotated.GetData()); 
	x.log_prob = x_rotated.log_prob; 
	x.SetWeight(x_rotated.GetWeight()); 
}

