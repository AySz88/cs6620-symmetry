#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm> // for std::max
#include <cstdlib> // for rand()
#include "pmvs_base/pmvs/patch.h"
#include "pmvs_base/numeric/vec4.h"

const int GRID_SIZE = 64;
const int SAMPLES = 1000;
const float THRESHOLD = 1000.0;
const int TOTAL_ISP_ITER = 20;
const int POINTS_PER_ITER = 1000;

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
					ReflectPlane p = getPlane(tIdx, pIdx, rIdx);
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

void main()
{
	patchVect patches;
	Vec4f avg;
	float maxR;
	readPatches(9, patches, avg, maxR); // avg has the average coord of patches BEFORE; it is subtracted out

	Vec3f minTPR(0, 0, 0), maxTPR(float(3.14159*0.5), float(2.0*3.14159), maxR);
	ReflectPlaneGrid grid(GRID_SIZE, minTPR, maxTPR);

	srand(1234); // DEBUG

	for(int sampNum = 0; sampNum < SAMPLES; sampNum++)
	{
		Cpatch p1 = patches[getRandom(0, patches.size())],
			p2 = patches[getRandom(0, patches.size())];

		Vec4f diff = p1.m_coord - p2.m_coord;
		float dSqrd = diff.norm2();
		float theta = acos(diff[2]);
		float sinTheta = sin(theta);
		float phi = acos(diff[0]/sinTheta);
		Vec3f normal(sinTheta*cos(phi), sinTheta*sin(phi), diff[2]);
		
		Vec4f intersectPt = float(0.5) * (p2.m_coord + p1.m_coord);
		Vec3f iPt3(intersectPt[0], intersectPt[1], intersectPt[2]);
		float r = iPt3 * normal;

		grid.getPlane(theta, phi, r).score += float(0.5)/(sinTheta*dSqrd);
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
		}
	}
	
	//ISP
	while(TOTAL_ISP_ITER)
	{
		for(std::vector<ReflectPlane>::iterator iter = symmetries.begin(); iter != symmetries.end(); iter++)
		{
			ReflectPlane candidateRefPlane = *iter;
			for(int i = 0; i < POINTS_PER_ITER; i++)
			{
				std::vector<Cpatch> patchArray;
				patchArray.push_back(patches[getRandom(0, patches.size())]);
			}
			
			for(int i = 0; i < POINTS_PER_ITER; i++)
			{
				std::vector<Cpatch> refPatchArray;
				
				//get the reflected point corresponding to each patch
				Vec3f refPoint;

				//find point in patches closest to refPoint;

				//add that point to refPatchArray


			}
			
			//find t p r for candidate plane that minimizes sum of distances between 
			//corresponding points in patchArray and refPatchArray

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