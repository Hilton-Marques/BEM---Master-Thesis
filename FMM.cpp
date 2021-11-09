#include<iostream>
#include "Solver.h"
#include "FMM.h"
#include "GaussQuadrature.h"
#include <algorithm> // set_difference
#include <set> // erase duplicates
extern const double M_PI;
extern const double Const;

FMM::FMM()
{
}
FMM::FMM(Solid* solid, int N, int NG, int Nc, int field) {


	m_nLMax = solid->m_nL;
	m_Nc = Nc;
	m_field = field;
	//m_elementsLevel = & solid->m_elementsLevel;
	m_nTotalVerts = solid->m_nPts;
	m_N = N;
	m_solver = new Solver(N, NG, m_nTotalVerts);
	m_solid = solid;
	m_nDb = m_solid->m_vertsBd.size(); // temperature prescribed
	m_Hdb = Matrix(m_nTotalVerts, m_nDb); 
	std::vector<Face*> elements;
	std::vector<Face*> elementsFather;
	elements = solid->m_elementsLevel[m_nLMax - 1];
	elementsFather = solid->m_elementsLevel[m_nLMax - 2];


	////Initialize  leaf elements
	m_tot = (N + 1) * (N + 1);

	for (Face* father : elementsFather)
	{
		father->setFMMSizes(m_tot);
		for (Face* element : father->m_childrenElement) 
		{
			element->m_points[0]->m_elements.push_back(element);
			element->m_points[1]->m_elements.push_back(element);
			element->m_points[2]->m_elements.push_back(element);
			element->m_adjacentFaces = element->getAdjacentElements();
			element->setNormal();
			element->m_isLeaf = true;
		}
	}

	double err;
	
	changeBoundaryCondition(solid->m_verts, elements);  // change boundary condition to potential fields


	m_b = computeBvector(err);


	//std::cout.precision(15);
	//std::cout << err << "\n";
	//std::cout << solid->m_nPts;
	//prepareForGMRES();
	//GMRES gmres(this);
	//m_x = gmres.Solver(m_b);
	//m_resid = gmres.m_resid;


	//Matrix A(4, 4, { 8, 7, 6, 1, 6, 7, 2, 5, 5, 9, 3, 2, 10, 9, 9, 10 });
	//Matrix b(4, 1, { 1,0,1,0 });
	//Matrix res = gmres.BiCGM(A,b);
	//for (int i = 0; i < 4; i++)
	//{
	//	std::cout << res[{i, 0}] << "\n";
	//}




		//checkWithConventionalBemH(solid, solid->m_elementsLevel.back());
		//checkWithConventionalBem(solid, solid->m_elementsLevel.back(), m_b);
	}



std::vector<Vertex*> FMM::getLeafNodesFromElements(std::vector<Face*> adjacentElements)
{
	std::vector<Vertex*> sourceVertss;
	std::vector<Vertex*> sourceVerts;
	if (adjacentElements.empty())
	{
		return sourceVerts;
	}	//sourceVerts.reserve(adjacentElements.size());

	for (Face* element : adjacentElements)
	{
		  std::vector<Face*> leafElements;
			element->getLeafElements(&leafElements);
			for (Face* leafElement : leafElements)
			{
				for (Vertex* pt : leafElement->m_points)
				{
					sourceVerts.push_back(pt);
				}
			}
	}
	// erase duplicates
	std::set<Vertex*> s(sourceVerts.begin(), sourceVerts.end());
	sourceVerts.assign(s.begin(), s.end());
	return sourceVerts;
}

std::vector<Face*> FMM::getFarElements(Face* element)
{
	std::vector<Face*> elements = m_solid->m_elementsLevel[element->m_nL - 1];
	std::vector<Face*> farElements;
	std::vector<Face*> adjElements = element->getAdjacentElements();
	adjElements.push_back(element);
	std::sort(elements.begin(), elements.end());
	std::sort(adjElements.begin(), adjElements.end());
	// got A \ B with poitns	
	std::set_difference(elements.begin(), elements.end(), adjElements.begin(), adjElements.end(), std::inserter(farElements, farElements.begin()));

	return farElements;
}

void FMM::checkWithConventionalBem(Solid* solid, std::vector<Face*> elements, Matrix &Ax)
{
	GaussQuad gauss(40);
	double* m_Ax;
	double* m_Bx;
	m_Ax = new double[solid->m_nPts];
	m_Bx = new double[solid->m_nPts];
	for (int i = 0; i < solid->m_nPts; i++)
	{
		m_Ax[i] = 0.0;
		m_Bx[i] = 0.0;
	}
	int count = 0;
#pragma omp parallel for
	for (int m = 0; m < elements.size(); m++)
	{
		Face* element = elements[m];
		//for (Face* element : elements) 

		Point* normal = &element->m_normal;
		double* jacobian = &element->m_jacobian;
		
		const double Const2 = Const  * (*jacobian);
		for (int i = 0; i < 3; i++)
		{
			const double q = element->m_q[i];
			for (Vertex* source : solid->m_verts)
			{
				double G = 0.0;
				// Integrate Gq
				for (int j = 0; j < gauss.m_totN; j++)
				{
					double* Xi = &gauss.m_gaussPointsXi[j];
					double* Eta = &gauss.m_gaussPointsEta[j];
					double* W = &gauss.m_gaussWeights[j];
					double shapeF[3] = { 1. - *Xi - *Eta  , *Xi, *Eta, };

					Point xksi = (shapeF[0] * element->m_points[0]->m_coord) + (shapeF[1] * element->m_points[1]->m_coord) + (shapeF[2] * element->m_points[2]->m_coord);
					Point r = xksi - source->m_coord;
					double rNorm = r.getNorm();
					G = G + (1 / rNorm) * (*W) * shapeF[i];
				}

				m_Ax[source->m_id] = m_Ax[source->m_id] + q * Const2 * G;
				m_Bx[source->m_id] = m_Bx[source->m_id] + q * Const2 * G;
			}
			count++;
		}
	}
	std::cout.precision(15);
	for (int i = 0; i < solid->m_nPts; i++)
	{
		std::cout << abs(m_Ax[i] - m_Bx[i]) << "\n";
	}
	delete[] m_Ax;
}

void FMM::checkWithConventionalBemH(Solid* solid, std::vector<Face*> elements)
{
	GaussQuad gauss(40);
	double* m_Ax;
	m_Ax = new double[solid->m_nPts];
	for (int i = 0; i < solid->m_nPts; i++)
	{
		m_Ax[i] = 0.0;
	}
	int count = 0;
#pragma omp parallel for
	for (int m = 0; m < elements.size(); m++)
	{
		Face* element = elements[m];
		//for (Face* element : elements) 

		Point* normal = &element->m_normal;
		double* jacobian = &element->m_jacobian;
		const int u = 1;
		const double Const2 = Const * u;

		for (int i = 0; i < 3; i++)
		{
			Vertex* pt = element->m_points[i];

			for (Vertex* source : solid->m_verts)
			{
				if (source->m_id == pt->m_id)
				{
					if (!source->markTemp)
					{
						// Calculate diagonal of H Matrix
						m_Ax[source->m_id] = m_Ax[source->m_id] + Const * source->calculateSolidAngle();
						source->markTemp = true;
					}
				}
				else {				
				double H = 0.0;
				// Integrate Gq
				for (int j = 0; j < gauss.m_totN; j++)
				{
					double* Xi = &gauss.m_gaussPointsXi[j];
					double* Eta = &gauss.m_gaussPointsEta[j];
					double* W = &gauss.m_gaussWeights[j];
					double shapeF[3] = { 1. - *Xi - *Eta  , *Xi, *Eta, };

					Point xksi = (shapeF[0] * element->m_points[0]->m_coord) + (shapeF[1] * element->m_points[1]->m_coord) + (shapeF[2] * element->m_points[2]->m_coord);
					Point r = xksi - source->m_coord;
					double rNorm = r.getNorm();

					double rNorm3 = rNorm * rNorm * rNorm;
					H = H + (-(r * (*normal)) / rNorm3) * (*W) * shapeF[i];

					
				}
				m_Ax[source->m_id] = m_Ax[source->m_id] + Const2 * H;
				}

				
			}
			count++;
		}
	}
	std::cout.precision(17);
	for (int i = 0; i < solid->m_nPts; i++)
	{
		std::cout << std::fixed << m_Ax[i] << "\n";
	}
	delete[] m_Ax;
}

Matrix FMM::computeBvector(double& erro)
{
	//tic();
	//Matrix Hd = this->Hd(m_solid, m_solid->m_elementsLevel[m_nLMax - 1]);
	//Matrix Gt = this->Gt(m_solid, m_solid->m_elementsLevel[m_nLMax - 1]);
	//toc();
	//Matrix bFMM = Matrix(m_nTotalVerts, 1);
	int Nc = m_Nc;
	int startLevel = m_nLMax - 1 - Nc;
	// Part 1

	std::vector<Face*> elements = m_solid->m_elementsLevel[startLevel];
	//tic();
	//conventional BEM

#pragma omp parallel for
		for (int i = 0; i < elements.size(); i++)
		{
			Face* element = elements[i];
			std::vector<Face*> adjacentElements = element->getAdjacentElements();
			std::vector<Vertex*> sourceVertexs = getLeafNodesFromElements(adjacentElements);
			std::vector<Face*> childrenElement;
			element->getLeafElements(&childrenElement);
			//element->setFMMSizes(m_tot);

			for (Face* child : childrenElement)
			{
				/*for (Vertex* source : sourceVertexs)
				{
					source->m_closeElements.push_back(child);
				}*/
				//m_solver->CalculateNearInt(sourceVertexs, child);
				m_solver->CalculateNearIntExact(sourceVertexs, child);
				// Integrate Expansion
				//m_solver->CalculateME(child);
				//child->m_collectVertexes.insert(child->m_collectVertexes.end(), sourceVertexs.begin(), sourceVertexs.end());
			}
		}
  //DEBUG
	//int nSource = 35;
	//std::vector<Face*> totalElements = m_solid->m_elementsLevel[m_nLMax - 1];
	//std::vector<Face*> closeElements = m_solid->m_verts[nSource]->m_closeElements;
	//std::vector<Face*> farElements;
	//std::sort(totalElements.begin(), totalElements.end());
	//std::sort(closeElements.begin(), closeElements.end());
	//std::set_difference(totalElements.begin(), totalElements.end(),
	//	closeElements.begin(), closeElements.end(), std::inserter(farElements, farElements.begin()));
	//double GtFmm = 0.0;
	//double HdFmm = 0.0;
	//double HdExact;
	//double GtExact;
	//for (Face* element : farElements)
	//{
	//	HdExact = m_solver->m_Hd[{nSource, 0}];
	//	GtExact = m_solver->m_Gt[{nSource, 0}];
 //   Point yc = element->m_fatherElement->getYc();
 //   std::vector<std::complex<double>> Sb;
 //   Point x = m_solid->m_verts[nSource]->m_coord - yc;
 //   Sb = m_solver->m_RTable.evaluateRecursiveTableS(x);
 //   double valueG = Const * m_solver->Dot(element->m_MEG, &Sb);
 //   double valueH = Const * m_solver->Dot(element->m_MEH, &Sb);
	//	m_solver->CalculateNearIntExact(std::vector<Vertex*> {m_solid->m_verts[nSource]}, element);
	//	double newHd = m_solver->m_Hd[{nSource, 0}] - HdExact;
	//	double newGt = m_solver->m_Gt[{nSource, 0}] - GtExact;
	//	double errHd = newHd - valueH;
	//	double errGt = newGt - valueG;
	//	GtFmm += valueG;
	//	HdFmm += valueH;
	//}
	//double err1 = m_solver->m_Gt[{nSource, 0}] - GtFmm;
	//double err2 = m_solver->m_Hd[{nSource, 0}] - HdFmm;


	// Part2
		if (Nc == 1)
		{
			for (int i = startLevel - 1; i > 0; i--)
			{
				elements = m_solid->m_elementsLevel[i];
				#pragma omp parallel for		
				for (int j = 0; j < elements.size(); j++)
				{
					Face* element = elements[j];

					//already evaluate FMM delivery downward pass
					std::vector<Face*> adjacentElementsMother = element->getAdjacentElements();				
					std::vector<Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsMother);
					std::sort(leafPoints.begin(), leafPoints.end());

					for (Face* child : element->m_childrenElement)
					{
						std::vector<Vertex*> sourceVertexs;
						std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
						std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
						std::vector<Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
						std::sort(adjacentNodes.begin(), adjacentNodes.end());
						// got A \ B with poitns	
						std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));		
						m_solver->CalculateFarInt(sourceVertexs, child);
					}
					element->setFMMSizes(m_tot);
					m_solver->CalculateMMT(element);					
				}
			}
    }
    else
    {
			for (int i = startLevel - 1; i > 0; i--)
			{
				elements = m_solid->m_elementsLevel[i];
				//#pragma omp parallel for		
				for (int j = 0; j < elements.size(); j++)
				{
					Face* element = elements[j];

					//already evaluate FMM delivery downward pass
					std::vector<Face*> adjacentElementsMother = element->getAdjacentElements();					
					std::vector<Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsMother);
					std::sort(leafPoints.begin(), leafPoints.end());

					for (Face* child : element->m_childrenElement)
					{
						std::vector<Vertex*> sourceVertexs;
						std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
						std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
						std::vector<Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
						std::sort(adjacentNodes.begin(), adjacentNodes.end());
						// got A \ B with poitns	
						std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));
						std::vector<Face*> granChildrenElement;
						child->getLevelElements(&granChildrenElement, i + Nc );
						for (Face* granChild : granChildrenElement)
						{
							for (Face* granGranChild : granChild->m_childrenElement)
							{
								m_solver->CalculateFarInt(sourceVertexs, granGranChild);
							}
							granChild->setFMMSizes(m_tot);
							m_solver->CalculateMMT(granChild);
						}
					}
				}
			}

    }

//Part 3
		if (startLevel != 0)
		{
			elements = m_solid->m_elementsLevel[0];
			std::vector<Vertex*> leafPoints = m_solid->m_verts;
			std::sort(leafPoints.begin(), leafPoints.end());

			//#pragma omp parallel for 	
			for (int i = 0; i < elements.size(); i++)
			{
				Face* element = elements[i];
				for (Face* child : element->m_childrenElement)
				{
					std::vector<Vertex*> sourceVertexs;
					std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
					std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
					std::vector<Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
					std::sort(adjacentNodes.begin(), adjacentNodes.end());
					// got A \ B with poitns	
					std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

					// Far integral
					//m_solver->CalculateFarInt(sourceVertexs, child);			
					std::vector<Face*> granChildrenElement;
					child->getLevelElements(&granChildrenElement, Nc + 1);
					for (Face* granChild : granChildrenElement)
					{
						//m_solver->CalculateFarInt(sourceVertexs, granChild);
						std::vector<Face*> grangranChildrenElement;
						granChild->getLeafElements(&grangranChildrenElement);
						Point yc = granChild->getYc();
						std::vector<std::complex<double>> Sb;
						for (Vertex* source : sourceVertexs)
						{
							for (Face* gran : grangranChildrenElement)
							{
								source->m_closeElements.push_back(gran);
							}
							Point x = source->m_coord - yc;
							Sb = m_solver->m_RTable.evaluateRecursiveTableS(x);
							double valueGt = Const * m_solver->Dot(granChild->m_MEG, &Sb);
							double valueHd = Const * m_solver->Dot(granChild->m_MEH, &Sb);
							m_solver->m_Gt[{source->m_id, 0}] += Const * m_solver->Dot(granChild->m_MEG, &Sb);
							m_solver->m_Hd[{source->m_id, 0}] += Const * m_solver->Dot(granChild->m_MEH, &Sb);
							
							//double valueGchild = 0.0;
							//double valueHchild = 0.0;
							//std::vector<Face*> leafElements;
							//child->getLeafElements(&leafElements);
							//for (Face* childEl : leafElements)
							//{
							//	x = source->m_coord - childEl->m_fatherElement->getYc();
							//	Sb = m_solver->m_RTable.evaluateRecursiveTableS(x);
							//	valueGchild += Const * m_solver->Dot(childEl->m_MEG, &Sb);
							//	valueHchild += Const * m_solver->Dot(childEl->m_MEH, &Sb);
							//}
							//double err1 = valueGt - valueGchild;
							//double err2 = valueHd - valueHchild;
						}
						granChild->reset();
					}
				}
			}
		}	
	 //toc();
	Matrix bFMM = m_solver->m_Gt - m_solver->m_Hd;
	//Matrix errHd = m_solver->m_Hd - Hd;
	//Matrix errGt = m_solver->m_Gt - Gt;
 // Matrix bBEM = Hd - Gt;
	//bFMM.show();
	//Matrix errHd = Hd - m_solver->m_Hd;
	//Matrix errGt = Gt - m_solver->m_Gt;
	//double value = err.norm() / Hd.norm();
	std::cout.precision(15);
	//std::cout << errHd.norm() << "\n";
	//std::cout << errGt.norm() << "\n";
	//std::cout << err.norm() << "\n";
	//std::cout << bFMM.norm() << "\n";
	
	double erroFMM = bFMM.norm() / m_solver->m_Hd.norm();
	//double erroBEM = bBEM.norm() / Hd.norm();
	//std::cout << m_solid->m_nPts << "\n";
	std::cout << "  "<< erroFMM << "\n";
	//std::cout << erroBEM << "\n";

	erro = bFMM.norm() / m_solver->m_Hd.norm();
	m_solver->reset();
	return bFMM;
}

Matrix FMM::matrixVectorMulti(Matrix &x)
{

	#pragma omp for
	for (int i = 0; i < m_nTotalVerts; i++)
	{
		Vertex* pt = m_solid->m_verts[i];
		if (pt->m_bd)
		{
			pt->m_u = x[{pt->m_id, 0}];
		}
	}

	// the same with elements
	std::vector<Face*> elements = m_solid->m_elementsLevel[m_nLMax - 2];
	//conventional BEM

#pragma omp parallel for 	

	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];
		std::vector<Face*> adjacentElements = element->getAdjacentElements();
		std::vector<Vertex*> sourceVertexs = getLeafNodesFromElements(adjacentElements);
		for (Face* child : element->m_childrenElement)
		{
			m_solver->CalculateNearIntMatrix(sourceVertexs, child);
			// Integrate Expansion
			//m_solver->CalculateME(child);

			for (int l = 0; l < 3; l++)
			{
				if (child->m_points[l]->m_bd)
				{
					for (int k = 0; k < m_tot; k++)
					{
						element->m_MEH[k] += child->m_points[l]->m_u * child->m_MEH_p[l][k];
					}
				}
				if (child->m_bd) // have to change to take care of three types of gradiente
				{
					for (int k = 0; k < m_tot; k++)
					{
						element->m_MEG[k] += child->m_q[l] * child->m_MEG_p[l][k];
					}
				}
			}

			//child->m_collectVertexes.insert(child->m_collectVertexes.end(), sourceVertexs.begin(), sourceVertexs.end());
			//child->reset();
		}

		// Accumulate moments

		//element->accumulate();

	}

	// Part2
	for (int i = m_nLMax - 3; i > 0; i--)
	{

		elements = m_solid->m_elementsLevel[i];
#pragma omp parallel for		
		for (int j = 0; j < elements.size(); j++)
		{
			Face* element = elements[j];

			//already evaluate FMM delivery downward pass
			std::vector<Face*> adjacentElementsMother = element->getAdjacentElements();
			std::vector<Face*> adjacentElementsChildren;
			adjacentElementsChildren.reserve(4 * adjacentElementsMother.size());
			for (Face* adj : adjacentElementsMother)
			{
				adjacentElementsChildren.insert(adjacentElementsChildren.end(), std::begin(adj->m_childrenElement), std::end(adj->m_childrenElement));
			}
			std::vector<Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsChildren);
			std::sort(leafPoints.begin(), leafPoints.end());


			for (Face* child : element->m_childrenElement)
			{
				std::vector<Vertex*> sourceVertexs;
				std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
				std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
				std::vector<Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
				std::sort(adjacentNodes.begin(), adjacentNodes.end());
				// got A \ B with poitns	
				std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

				//far integral
				m_solver->CalculateFarInt(sourceVertexs, child);

			}
			element->setFMMSizes(m_tot);
			std::vector<Face*> childrenElement = element->getChildrenElement();
			m_solver->CalculateMMT(element);
		}
	}
	// Part3 ///////////////////
	elements = m_solid->m_elementsLevel[0];
//#pragma omp parallel for 	
	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];

		//already evaluate FMM delivery downward pass
		std::vector<Face*> adjacentElementsMother = element->m_adjacentFaces;
		std::vector<Face*> adjacentElementsChildren;
		adjacentElementsChildren.reserve(4 * adjacentElementsMother.size());
		for (Face* adj : adjacentElementsMother)
		{
			adjacentElementsChildren.insert(adjacentElementsChildren.end(), std::begin(adj->m_childrenElement), std::end(adj->m_childrenElement));
		}
		std::vector<Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsChildren);
		std::sort(leafPoints.begin(), leafPoints.end());

		auto farElements = getFarElements(element);
		auto farNodes = getLeafNodesFromElements(farElements);

		for (Face* child : element->m_childrenElement)
		{
			std::vector<Vertex*> sourceVertexs;
			std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
			std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
			std::vector<Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
			std::sort(adjacentNodes.begin(), adjacentNodes.end());
			// got A \ B with poitns	
			std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

			// accumulate nodes from far elements
			if (!farNodes.empty())
			{
				sourceVertexs.insert(sourceVertexs.end(), farNodes.begin(), farNodes.end());
			}

			// Far integral
			//m_solver->CalculateFarInt(sourceVertexs, child );

			Point yc = child->getYc();
			std::vector<std::complex<double>> Sb;
			for (Vertex* source : sourceVertexs)
			{

				Point x = source->m_coord - yc;
				Sb = m_solver->m_RTable.evaluateRecursiveTableS(x);
				m_solver->m_Gt[{source->m_id, 0}] += Const * m_solver->Dot(child->m_MEG, &Sb);
				m_solver->m_Hd[{source->m_id, 0}] += Const * m_solver->Dot(child->m_MEH, &Sb);

			}
			child->reset();
		}
	}
	for (int i = 0; i < m_nDb; i++)
	{
		Vertex* pt = m_solid->m_vertsBd[i];
		m_solver->m_Hd[{pt->m_id, 0}] = x[{pt->m_id, 0}];
	}
	Matrix res = (m_solver->m_Hd - m_solver->m_Gt);
	m_solver->reset();
	return res;
}

void FMM::prepareForGMRES()
{
	m_solver->m_type = SolverType::GMRES;
	// inverte boundary conditions
	for (Vertex* pt : m_solid->m_verts)
	{
		if (pt->m_bd)
		{
			pt->m_bd = false;
		}
		else
		{
			pt->m_bd = true;
			pt->m_u = 1.0; // calculate ME temporary
		}
		pt->markTemp = false;
	}

	std::vector<Face*> elements = m_solid->m_elementsLevel[m_nLMax - 1];

	for (Face* element : elements)
	{
		//initialize vectors 
		element->set_FMM_ME_sizes();
		if (element->m_bd)
		{
			element->m_bd = false;
		}
		else
		{
			for (int i = 0; i < 3; i++)
			{
				element->m_bd = true;
				element->m_q[i] = 1.0; // calculate ME temporary

			}

		}
	}

	elements = m_solid->m_elementsLevel[m_nLMax - 2];

#pragma omp parallel for		
	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];
		std::vector<Face*> adjacentElements = element->getAdjacentElements();
		std::vector<Vertex*> sourceVertexs = getLeafNodesFromElements(adjacentElements);
		for (Face* child : element->m_childrenElement)
		{	
			m_solver->prepareElementMatrix(sourceVertexs, child);		
			m_solver->CalculateMEForGMRES(child);
		}
	}


}

void FMM::showB()
{
	std::cout.precision(15);
	for (int i = 0; i < m_nTotalVerts; i++)
	{
		std::cout << m_b[{i, 0}] << "\n";
	}
}

void FMM::showX()
{
	std::cout.precision(15);
	for (int i = 0; i < m_nTotalVerts; i++)
	{
		std::cout << m_x[{i, 0}] << "\n";
	}
}

double FMM::trivialField(double x, double y, double z)
{
	switch (m_field)
	{
	default:
		break;
	case 0:
		//return (x + y);
		return (x + y + z);
	case 1:
		//return ((x * x) - (z * z));
		return pow(x, 2) + pow(y, 2) - 2 * pow(z, 2);
	case 2:
		return (x * y * z);
	}

	//return x+y;	
	//return ((x * x) - (z * z));
	//return (pow(x, 3) * z - 3 * x * pow(y, 2) * z);
	//-(Power(x,2)*z) + Power(y,2)*z
}

double FMM::trivialGrad(double x, double y, double z, Point& n)
{
	switch (m_field)
	{
	default:
		break;
	case 0:
	{
		//Point grad_linear(1, 1, 0);
		Point grad_linear(1, 1, 1);
		return grad_linear * n;
	}
	case 1:
	{
		//Point grad_quad(2 * x, 0, -2 * z);
		Point grad_quad(2 * x, 2 * y, -4 * z);
		return grad_quad * n;
	}
	case 2:
	{
		Point grad_cubic(y * z, x * z, x * y);
		return grad_cubic * n;
	}
	}
	//Point grad(1, 1, 1);
	//Point grad(1, 1, 0);
	//Point grad(2 * x, 0, -2 * z);
	//Point grad(3 * pow(x, 2) * z - 3 * pow(y, 2) * z, -6 * x * y * z, pow(x, 3) - 3 * x * pow(y, 2));
	//return grad * n;
	//
}

void FMM::changeBoundaryCondition(std::vector<Vertex*> points, std::vector<Face*> elements)
{
	for (Vertex* pt : points)
	{
		pt->m_u = trivialField(pt->m_coord.m_x, pt->m_coord.m_y, pt->m_coord.m_z);
		//pt->m_u = 1.0; // diagonal
		pt->m_bd = true;
	}
	for (Face* element : elements)
	{
		element->m_bd = true;
		Point n = element->getNormalUnitary();
		double value = 0.0;
		for (int i = 0; i < 3; i++)
		{
			Point yc = element->m_points[i]->m_coord; 
			value += trivialGrad(yc.m_x, yc.m_y, yc.m_z, n);
		}
		value = value / 3;
		element->m_q[0] = value;
		element->m_q[1] = value;
		element->m_q[2] = value;
		//element->m_q[0] = 0;
		//element->m_q[1] = 0;
		//element->m_q[2] = 0;
	}
}

void FMM::toc()
{
	std::cout << " "
		<< ((double)(clock() - m_tictoc_stack.top())) / CLOCKS_PER_SEC
		<< std::endl;
	m_tictoc_stack.pop();
}





FMM::~FMM()
{

}

Matrix FMM::Hd(Solid* solid, std::vector<Face*> elements)
{
	Matrix Hd(solid->m_nPts, 1);
	int count = 0;
	#pragma omp parallel for
	for (int m = 0; m < elements.size(); m++)
	{
		Face* element = elements[m];

		for (Vertex* source : solid->m_verts)
		{
			Matrix integral = element->exactIntH(source);
			for (int i = 0; i < 3; i++)
			{
				Vertex* pt = element->m_points[i];
				double Const2 = Const * pt->m_u;

				if (source->m_id == pt->m_id)
				{
					if (!source->markTemp)
					{
						// Calculate diagonal of H Matrix
						source->markTemp = true;
						#pragma omp atomic
						Hd[{source->m_id, 0}] += Const2 * source->calculateSolidAngle();
					}
				}
				else
				{
					#pragma omp atomic
					Hd[{source->m_id, 0}] += Const2 * integral[{i, 0}];
				}
			}
		}
	}

	for (Vertex* source : solid->m_verts)
	{
		source->markTemp = false;
	}
	return Hd;

}


Matrix FMM::Gt(Solid* solid, std::vector<Face*> elements)
{
	Matrix Gt(solid->m_nPts, 1);
	int count = 0;
#pragma omp parallel for
	for (int m = 0; m < elements.size(); m++)
	{
		Face* element = elements[m];
		double* jacobian = &element->m_jacobian;
		const double Const2 = Const  * (*jacobian);
		for (Vertex* source : solid->m_verts)
		{
			Matrix integral = element->exactIntG(source);
			for (int i = 0; i < 3; i++)
			{
				double q = element->m_q[i];
				#pragma omp atomic
				Gt[{source->m_id, 0}] += q * Const2 * integral[{i, 0}];
			}
		}
	}
	return Gt;

}