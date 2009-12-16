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
const float THRESHOLD = 0.0001;
const int TOTAL_ISP_ITER = 20;
const int POINTS_PER_ITER = 1000;
const float SIGMA = 1;
const float SIGMASQRD = SIGMA*SIGMA;
const float PI = 3.1415926535;

typedef TVec2<float> Vec2f;
typedef TVec3<float> Vec3f;
typedef TVec4<float> Vec4f;

typedef Patch::Cpatch Cpatch;
typedef std::vector<Cpatch> patchVect;

void readPsetFile(patchVect& patches, Vec4f& avg, float& maxR);

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
	if (theta > 0.5*PI)
	{
		diff *= -1;
		theta = acos(diff[2]);
		needNegZ = true;
	}
	sinTheta = sin(theta);

	// FIXME something is wrong here?
	phi = acos(std::max<float>(std::min<float>(diff[0]/sinTheta, 1.0), -1.0));
	if (diff[1] < 0) phi = 2*PI - phi;
	//phi = acos(std::max<float>(std::min<float>(diff[0]/sinTheta, 1.0), 0.0)); // guard against rounding errors
	Vec3f normal(sinTheta*cos(phi), sinTheta*sin(phi), diff[2]);
	
	// DEBUG
	/*
	float normalError = (normal - Vec3f(diff[0], diff[1], diff[2])).norm();
	if (normalError > 1e-3)
	{
		std::cout << "Error in normal was" << normalError << std::endl;
		std::cout << "  Normal:" << normal[0] << " " << normal[1] << " " << normal[2] << " " << normal[3] << std::endl;
		std::cout << "  Diff:" << diff[0] << " " << diff[1] << " " << diff[2] << " " << diff[3] << std::endl;
	}*/
	
	Vec4f intersectPt = float(0.5) * (p2.m_coord + p1.m_coord);
	Vec3f iPt3(intersectPt[0], intersectPt[1], intersectPt[2]);
	r = iPt3 * normal * (needNegZ ? -1.0 : 1.0);

	if (theta > 0.5*PI || theta < 0 || phi < 0 || phi > 2*PI)
	{
		/*
		r *= -1;
		theta = PI - theta;
		phi = PI + phi;
		if (phi > 2.0*PI) phi -= 2.0*PI;
		*/
		std::cout << "Bad reflection calculated! " << p1.m_coord - p2.m_coord << " --> "
			<< theta << " " << phi << " " << r << std::endl;
	}
}

void writePlyFiles(patchVect patches, std::vector<ReflectPlane> symmetries)
{
	int i = 0;
	
	
	for(std::vector<ReflectPlane>::iterator iter = symmetries.begin(); iter != symmetries.end(); iter++)
	{	
		std::stringstream filename;
		filename<<i<<".ply";
		i++;
		std::ofstream data;
		data.open(filename.str().c_str());
		ReflectPlane curPlane = *iter;
		data<<"ply\nformat ascii 1.0\nelement vertex "<<patches.size()<<"\nproperty float x\nproperty float y\n"
				"property float z\nproperty uchar diffuse_red\nproperty uchar diffuse_green\n"
				"property uchar diffuse_blue\nend_header"<<std::endl;
		for(patchVect::iterator iter1 = patches.begin(); iter1 != patches.end(); iter1++)
		{
			Cpatch p = *iter1;
			float  sinTheta = sin(curPlane.theta);
			Vec4f normal(sinTheta*cos(curPlane.phi), sinTheta*sin(curPlane.phi), cos(curPlane.theta), 0.0);
			float dotp = p.m_coord * normal;
			std::stringstream line;
			if(dotp < 0.0) // FIXME if(dotp < curPlane.r) ?
			{
				data<<p.m_coord[0]<<" "<<p.m_coord[1]<<" "<<p.m_coord[2]<<" "<<0<<" "<<255<<" "<<0<<std::endl;
				
			}
			else
			{
				data<<p.m_coord[0]<<" "<<p.m_coord[1]<<" "<<p.m_coord[2]<<" "<<255<<" "<<0<<" "<<0<<std::endl;
			}
			
		}
		data.close();
	}
}


void main()
{
	patchVect patches;
	Vec4f avg;
	float maxR;
	//readPatches(9, patches, avg, maxR); // avg has the average coord of patches BEFORE; it is subtracted out
	readPsetFile(patches, avg, maxR);
	/*
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
	*/
	Vec3f minTPR(0, 0, -maxR), maxTPR(float(PI*0.5), float(2.0*PI), maxR);
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
	
	writePlyFiles(patches, symmetries);

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

void readPsetFile(patchVect& patches, Vec4f& avg, float& maxR)
{
		std::stringstream filename;
		filename << "mcgraw.pset";

		std::cout << "Starting to read image " << filename.str() << std::endl;

		std::ifstream data;
		data.open(filename.str().c_str());
		int points;
		std::string line;
		
			
		while(1)
		{
			getline(data,line);
			if (line[0] != 'E')
			{
				int i =0, j = 0;
				Cpatch p;
				char * coordStrX = new char[20];
				char * coordStrY = new char[20];
				char * coordStrZ = new char[20];
				while(line[i] != ' ')
				{
					coordStrX[j++] = line[i++];
				}
				coordStrX[j] = '\0';
				j=0;
				i++;
				while(line[i] != ' ')
				{
					coordStrY[j++] = line[i++];
				}
				coordStrY[j] = '\0';
				j =0;
				i++;
				while(line[i] != ' ')
				{
					coordStrZ[j++] = line[i++];
				}
				coordStrZ[j] = '\0';
				p.m_coord[0] = atof(coordStrX);
				p.m_coord[1] = atof(coordStrY);
				p.m_coord[2] = atof(coordStrZ);
				p.m_coord[3] = 1.0;
				patches.push_back(p);
				avg += p.m_coord;

			}
			else
				break;
		}
		std::cout<< "Done reading data from file"<<std::endl;

	avg /= avg[3];
	maxR = 0.0f;
	for(patchVect::iterator iter = patches.begin(); iter != patches.end(); iter++)
	{
		iter->m_coord -= avg;
		maxR = std::max(maxR, iter->m_coord.norm2());
	}
	maxR = sqrt(maxR);
}

volatile int getRandom(int min, int max)
{
	return (rand() % (max - min)) + min;
}

volatile float getRandomF(float min, float max)
{
	return float(rand())/float(RAND_MAX) * (max - min) + min;
}