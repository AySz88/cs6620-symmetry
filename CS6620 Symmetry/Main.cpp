#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm> // for std::max
#include <cstdlib> // for rand()
#include <limits>
#include "pmvs_base/pmvs/patch.h"
#include "pmvs_base/numeric/vec4.h"

const int GRID_SIZE = 64;
const int SAMPLES = 1000;
const float THRESHOLD = 10.0;
const int TOTAL_ISP_ITER = 20;
const int POINTS_PER_ITER = 1000;
const float SIGMA = 1;
const float SIGMASQRD = SIGMA*SIGMA;

typedef TVec2<float> Vec2f;
typedef TVec3<float> Vec3f;
typedef TVec4<float> Vec4f;

typedef Patch::Cpatch Cpatch;
typedef std::vector<Cpatch> patchVect;

struct ReflectPlane
{
public:
	float theta, phi, r; // lower corner (min_theta, min_phi, etc)
	float score;
};

class ReflectPlaneGrid
{
	ReflectPlane* grid;
	int size; // per dimension
	Vec3f min, max, diff; // theta, phi, and r
public:
	ReflectPlaneGrid(int gridSize, Vec3f minTPR, Vec3f maxTPR)
	{
		size = gridSize;
		grid = new ReflectPlane[size*size*size];
		min = minTPR;
		max = maxTPR;
		diff = (max - min)/size;

		float theta = min[0];
		for(int tIdx = 0; tIdx < size; tIdx++)
		{
			float phi = min[1];
			for(int pIdx = 0; pIdx < size; pIdx++)
			{
				float r = min[2];
				for (int rIdx = 0; rIdx < size; rIdx++)
				{
					ReflectPlane& p = getPlane(tIdx, pIdx, rIdx);
					p.theta = theta;
					p.phi = phi;
					p.r = r;
					p.score = 0;

					r += diff[2];
				}

				phi += diff[1];
			}

			theta += diff[0];
		}
	}

	~ReflectPlaneGrid()
	{
		delete[] grid;
	}

	ReflectPlane& getPlane(int t, int p, int r)
	{
		return grid[((t*size)+p)*size+r];
	}

	ReflectPlane& getPlane(float t, float p, float r)
	{
		return getPlane(int((t-min[0])/diff[0]), int((p-min[1])/diff[1]), int((r-min[2])/diff[2]));
	}

	bool isLocalMax(int t, int p, int r) // score
	{
		float thisScore = getPlane(t,p,r).score;
		return !(
			existsAndHigherScore(thisScore, t,p,r-1)
			|| existsAndHigherScore(thisScore, t,p,r+1)
			|| existsAndHigherScore(thisScore, t,p-1,r-1)
			|| existsAndHigherScore(thisScore, t,p-1,r)
			|| existsAndHigherScore(thisScore, t,p-1,r+1)
			|| existsAndHigherScore(thisScore, t,p+1,r-1)
			|| existsAndHigherScore(thisScore, t,p+1,r)
			|| existsAndHigherScore(thisScore, t,p+1,r+1)
			|| existsAndHigherScore(thisScore, t-1,p-1,r-1)
			|| existsAndHigherScore(thisScore, t-1,p-1,r)
			|| existsAndHigherScore(thisScore, t-1,p-1,r+1)
			|| existsAndHigherScore(thisScore, t-1,p,r-1)
			|| existsAndHigherScore(thisScore, t-1,p,r)
			|| existsAndHigherScore(thisScore, t-1,p,r+1)
			|| existsAndHigherScore(thisScore, t-1,p+1,r-1)
			|| existsAndHigherScore(thisScore, t-1,p+1,r)
			|| existsAndHigherScore(thisScore, t-1,p+1,r+1)
			|| existsAndHigherScore(thisScore, t+1,p-1,r-1)
			|| existsAndHigherScore(thisScore, t+1,p-1,r)
			|| existsAndHigherScore(thisScore, t+1,p-1,r+1)
			|| existsAndHigherScore(thisScore, t+1,p,r-1)
			|| existsAndHigherScore(thisScore, t+1,p,r)
			|| existsAndHigherScore(thisScore, t+1,p,r+1)
			|| existsAndHigherScore(thisScore, t+1,p+1,r-1)
			|| existsAndHigherScore(thisScore, t+1,p+1,r)
			|| existsAndHigherScore(thisScore, t+1,p+1,r+1)
			);
	}

private:
	bool existsAndHigherScore(float score, int t, int p, int r)
	{
		return
			t >= 0
			&& t < size
			&& p >= 0
			&& p < size
			&& r >= 0
			&& r < size
			&& getPlane(t,p,r).score > score;
	}
};

// helper
void readPatches(int files, patchVect& patches, Vec4f& avg, float& maxR);
volatile int getRandom(int min, int max);
volatile float getRandomF(float min, float max);
void getPlaneReflecting(Cpatch& p1, Cpatch& p2, float& theta, float& phi, float& r, float& sinTheta, float& dSqrd)
{
	Vec4f diff = p1.m_coord - p2.m_coord;
	dSqrd = diff.norm2();
	diff[3] = 0.0;
	diff.unitize();
	theta = acos(diff[2]);
	bool needNegZ = false;
	if (theta > 0.5*3.14159)
	{
		diff *= -1;
		theta = acos(diff[2]);
		needNegZ = true;
	}
	sinTheta = sin(theta);
	phi = acos(std::max<float>(std::min<float>(diff[0]/sinTheta, 1.0), 0.0)); // guard against rounding errors
	Vec3f normal(sinTheta*cos(phi), sinTheta*sin(phi), diff[2]);
	
	Vec4f intersectPt = float(0.5) * (p2.m_coord + p1.m_coord);
	Vec3f iPt3(intersectPt[0], intersectPt[1], intersectPt[2]);
	r = iPt3 * normal * (needNegZ ? -1.0 : 1.0);

	if (theta > 0.5*3.14159 || theta < 0)
	{
		/*
		r *= -1;
		theta = 3.14159 - theta;
		phi = 3.14159 + phi;
		if (phi > 2.0*3.14159) phi -= 2.0*3.14159;
		*/
		std::cout << "Bad reflection calculated! " << p1.m_coord - p2.m_coord << " --> "
			<< theta << " " << phi << " " << r << std::endl;
	}
}

void main()
{
	patchVect patches;
	Vec4f avg;
	float maxR;
	//readPatches(9, patches, avg, maxR); // avg has the average coord of patches BEFORE; it is subtracted out

	for(int i = 0; i < 10000; i++)
	{
		Cpatch p;
		p.m_coord[0] = getRandomF(-1.0, 1.0);
		p.m_coord[1] = getRandomF(-1.0, 1.0);
		p.m_coord[2] = getRandomF(-1.0, 1.0);
		p.m_coord[3] = 0.0;
		p.m_coord.unitize();
		p.m_coord *= getRandomF(0.9, 1.1);

		patches.push_back(p);
	}
	maxR = 1.1;
	avg[0] = 0.0;
	avg[1] = 0.0;
	avg[2] = 0.0;
	avg[3] = 1.0;

	Vec3f minTPR(0, 0, -maxR), maxTPR(float(3.14159*0.5), float(2.0*3.14159), maxR);
	ReflectPlaneGrid grid(GRID_SIZE, minTPR, maxTPR);

	srand(1234); // DEBUG

	for(int sampNum = 0; sampNum < SAMPLES; sampNum++)
	{
		Cpatch p1 = patches[getRandom(0, patches.size())],
			p2 = patches[getRandom(0, patches.size())];

		float theta, phi, r, sinTheta, dSqrd;
		getPlaneReflecting(p1, p2, theta, phi, r, sinTheta, dSqrd);

		if (sinTheta*dSqrd < 0 || sinTheta*dSqrd != sinTheta*dSqrd)
		{
			std::cout << "Negative score: 0.5/(" << sinTheta << "*" << dSqrd << " = " << float(0.5)/(sinTheta*dSqrd)
				<< " at " << theta << " " << phi << " " << r << std::endl;
		}

		grid.getPlane(theta, phi, r).score += float(0.5)/(sinTheta*dSqrd);

		/*
		std::cout << "score: " << grid.getPlane(theta, phi, r).score << std::endl;
		std::cout << "  added 0.5/(" << sinTheta << "*" << dSqrd << " = " << float(0.5)/(sinTheta*dSqrd)
			<< " at " << theta << " " << phi << " " << r << std::endl;
		*/
	}

	ReflectPlane* maxSoFar = &grid.getPlane(0,0,0);
	std::vector<ReflectPlane> localMaxes;
	for(int x = 0; x < GRID_SIZE; x++)
	{
		for(int y = 0; y < GRID_SIZE; y++)
		{
			for(int z = 0; z < GRID_SIZE; z++)
			{
				ReflectPlane& p = grid.getPlane(x,y,z);
				//if (p.score > 0) std::cout << "Score = " << p.score;
				if (p.score > THRESHOLD * (1.0 - p.r / maxR)
					&& p.score > 0.1*maxSoFar->score
					&& grid.isLocalMax(x,y,z)
					)
				{
					localMaxes.push_back(p);
					if (p.score > maxSoFar->score)
					{
						maxSoFar = &p;
					}
				}
			}
		}
	}

	std::vector<ReflectPlane> symmetries;
	// keep only large relative to max found
	for(std::vector<ReflectPlane>::iterator iter = localMaxes.begin(); iter != localMaxes.end(); iter++)
	{
		if (iter->score > 0.1*maxSoFar->score)
		{
			symmetries.push_back(*iter);
			//std::cout << *iter;
			std::cout << "Added plane " << iter->theta << " " << iter->phi << " " << iter->r << std::endl;
		}
	}
	
	//ISP
	while(TOTAL_ISP_ITER)
	{
		std::vector<Cpatch> patchArray;
		for(int i = 0; i < POINTS_PER_ITER; i++)
		{
			patchArray.push_back(patches[getRandom(0, patches.size())]);
		}

		for(std::vector<ReflectPlane>::iterator iter = symmetries.begin(); iter != symmetries.end(); iter++)
		{
			ReflectPlane candidateRefPlane = *iter;
			
			std::vector<Cpatch> refPatchArray;
			for(int i = 0; i < POINTS_PER_ITER; i++)
			{
				Vec4f curPoint = patchArray[i].m_coord;

				//get the reflected point corresponding to each patch
				float sinTheta = sin(iter->theta);
				float phi = iter->phi;
				Vec4f normal(sinTheta*cos(iter->phi), sinTheta*sin(iter->phi), cos(iter->theta), 0.0);

				float dist = sqrt(normal * curPoint);
				Vec4f projection = dist * normal;
				Vec4f perp = curPoint - projection;
				float reflDist = 2*iter->r - dist;

				Vec4f refPoint = reflDist * normal + perp;

				//find point in patches closest to refPoint
				// FIXME use octree?
				Cpatch& minSoFar = patches[0];
				float minDistSqrdSoFar = (patches[0].m_coord - refPoint).norm2();
				for (std::vector<Cpatch>::iterator iter2 = patches.begin(); iter2 != patches.end(); iter2++)
				{
					float thisDist = (iter2->m_coord - refPoint).norm2();
					if (thisDist < minDistSqrdSoFar)
					{
						minDistSqrdSoFar = thisDist;
						minSoFar = *iter2;
					}
				}

				//add that point to refPatchArray
				refPatchArray.push_back(minSoFar);
			}
			
			//find t p r for candidate plane that minimizes sum of distances between 
			//corresponding points in patchArray and refPatchArray
			
			float tTot = 0, pTot = 0, rTot = 0, weight = 0;
			for (int i = 0; i < POINTS_PER_ITER; i++)
			{
				if ((patchArray[i].m_coord - refPatchArray[i].m_coord).norm2() > 1e-5)
				{
					float t, p, r, sinT, dSqrd;
					getPlaneReflecting(patchArray[i], refPatchArray[i], t, p, r, sinT, dSqrd);
					tTot += t;
					pTot += p;
					rTot += r;
					weight += exp(-dSqrd/SIGMASQRD);

					if (p!=p || r!=r)
					{
						std::cout << "Weight " << exp(-dSqrd/SIGMASQRD)
							<< "  -  Plane " << iter->theta << " " << iter->phi << " " << iter->r
							<< " tugged toward " << t << " " << p << " " << r << std::endl;
						
						getPlaneReflecting(patchArray[i], refPatchArray[i], t, p, r, sinT, dSqrd);
					}
				}
			}
			tTot /= weight;
			pTot /= weight;
			rTot /= weight;

			std::cout << "Plane " << iter->theta << " " << iter->phi << " " << iter->r
				<< " updated to " << tTot << " " << pTot << " " << rTot << std::endl;

			// update to mean values
			iter->theta = tTot;
			iter->phi = pTot;
			iter->r = rTot;
		}

	}//ISP

	std::cout << "Done! type some char to exit..." << std::endl;
	char c;
	std::cin >> c;
}

void readPatches(int files, patchVect& patches, Vec4f& avg, float& maxR)
{
	for(int i = 0; i < files; i++)
	{
		// file names expected to be 0000.patc0, 0001.patc0, etc.
		std::stringstream filename;
		filename << std::setw(4) << std::setfill('0') << i << ".patc0";

		std::cout << "Starting to read image " << filename.str() << std::endl;

		std::ifstream data;
		data.open(filename.str().c_str());

		if (data.is_open())
		{
			std::string line;

			// value of "A" from PMVS-1's match.exe (ignored)
			getline(data, line);

			// number of points in file
			int points;
			data >> points;
			
			std::cout << "Points in image " << i << ": " << points << std::endl;

			//int count = 0;
			while (data.good())//&& count < points)
			{
				// skip lines until find "PATCH0"
				//std::streampos start;
				do
				{
					//start = data.tellg();
					getline(data, line);

					//std::cout << "line is " << line << std::endl;
				}
				while (data.good() && !(line == "PATCH0"));

				if (data.good())
				{
					//data.seekg(start); // for some reason this doesn't work, so just telling patch.>> to not look for the header

					// read in
					Patch::Cpatch p;
					data >> p;
					patches.push_back(p);
					avg += p.m_coord;
					
					//std::cout << p << std::endl;
				}
				//count++;
			}
			data.close();
		}
		else
		{
			std::cout << "Couldn't read file!" << std::endl;
		}
	} // end loop over files

	avg /= avg[3];
	maxR = 0.0f;
	for(patchVect::iterator iter = patches.begin(); iter != patches.end(); iter++)
	{
		iter->m_coord -= avg;
		maxR = std::max(maxR, iter->m_coord.norm2());
	}
}

volatile int getRandom(int min, int max)
{
	return (rand() % (max - min)) + min;
}

volatile float getRandomF(float min, float max)
{
	return float(rand())/float(RAND_MAX) * (max - min) + min;
}