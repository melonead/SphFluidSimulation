#include "Hashtable.h"
#include <cstdint>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <cassert>

HashTable::HashTable()
{

	// settings

	settings = SettingsSingleton::instance();
	settings->test = 3.0f;

	numberOfCells = (settings->containerWidth / settings->cellSize) * (settings->containerHeight / settings->cellSize);

	// Add an extra cell so that we don't index outside the array while computing
	// partial sums.
	grid = std::vector<unsigned int>(numberOfCells + 1, 0);

	query.neighborBucket.reserve(72);
	query.neighborBucket.resize(72);

	mouseQuery.neighborBucket.reserve(144);
	mouseQuery.neighborBucket.resize(144);
}

void HashTable::generateParticlesIDs(unsigned int numParticles)
{
	particleIDs = std::vector<unsigned int>(numParticles, 0);
	sortedParticleIDs = std::vector<unsigned int>(numParticles, 0);
	particleCounts = std::vector<unsigned int>(numberOfCells, 0);

	// create particles IDs
	for (unsigned int i = 0; i < numParticles; i++)
	{
		particleIDs[i] = i;
	}
}



NeighborQuery HashTable::getNeighborIDs(glm::vec2& position)
{
	// REVISIT: Perhaps make this a member variable?
	unsigned int neighborCellsCount = 9;

	// REVISIT:(temporary)
	// NeighborQuery query;
	// query.neighborBucket.resize(72);

	// We are creating this here insteading of using a single one because
	// of multithreading. 
	NeighborQuery query;
	query.size = 0;
	query.neighborBucket.resize(72);

	glm::vec2 currentCellPos = getCellPosition(position);
	glm::vec2 neighborCellPos;
	glm::vec2 pos;
	unsigned int key;
	for (unsigned int i = 0; i < neighborCellsCount; i++) {
		// REVISIT: Do the cell displacements need to be scaled by
		// the cell size?

		pos.x = position.x + cellDisplacements[i][0] * settings->cellSize;
		pos.y = position.y + cellDisplacements[i][1] * settings->cellSize;
		neighborCellPos = getCellPosition(pos);
		key = tableHash(neighborCellPos);
		// Go to table and find where to start looking
		unsigned int start = grid[key];
		// Go to particleCount and find how much to advance start at start
		unsigned int size = particleCounts[key];

		for (unsigned int j = 0; j < size; j++)
		{
			
			unsigned int id = sortedParticleIDs[start + j];
			 if (query.size >= query.neighborBucket.size())
			 {
			}
			query.neighborBucket.resize(query.size + 10);

			query.neighborBucket[query.size] = id;
			query.size += 1;
		}
	}

	return query;
}



NeighborQuery& HashTable::getNeighborIDsForMouse(glm::vec2& position, float radius)
{
	float startX = 0.0f;
	float startY = 0.0f;

	unsigned int boundingWidth = 9;
	mouseQuery.size = 0;

	glm::vec2 pos;
	glm::vec2 cellPos;
	// std::cout << "---------------------------------------------------------------------" << std::endl;
	for (unsigned int y = 0; y < boundingWidth; y++)
	{
		for (unsigned int x = 0; x < boundingWidth; x++)
		{

			pos.x = position.x + x * settings->cellSize;
			pos.y = position.y + y * settings->cellSize;

			//std::cout << "x: " << pos.x << " y: " << pos.y << std::endl;
			
			cellPos = getCellPosition(pos);

			unsigned int key = computeKey(cellPos);


			// Go to table and find where to start looking
			unsigned int start = grid[key];
			// Go to particleCount and find how much to advance start at start
			unsigned int size = particleCounts[key];

			for (unsigned int k = 0; k < size; k++)
			{
				unsigned int id = sortedParticleIDs[start + k];

				if (mouseQuery.size >= mouseQuery.neighborBucket.size())
				{
					mouseQuery.neighborBucket.resize(mouseQuery.size + 10);
				}

				mouseQuery.neighborBucket[mouseQuery.size] = id;
				mouseQuery.size += 1;
			}


		}
	}
	return mouseQuery;
}





void HashTable::createTable(std::vector<glm::vec2>& positions)
{
	grid = std::vector<unsigned int>(numberOfCells + 1, 0);
	particleCounts = std::vector<unsigned int>(numberOfCells, 0);
	unsigned int numParticles = positions.size();
	// Compute cell particle count
	for(unsigned int i = 0; i < numParticles; i++)
	{
		glm::vec2& position = positions[i];
		unsigned int key = computeKey(position);

		grid[key] += 1;
		// This apparent duplication is necessary because cellParticlesCount is going to change.
		particleCounts[key] += 1;
	}
	// compute partial sums
	int k = 0;
	for(unsigned int i = 0; i < numberOfCells; i++)
	{
		grid[i + 1]  += grid[i];
	}
	// Sort the particle ID, particles in the same cell will be next to each other
	
	for(unsigned int i = 0; i < numParticles; i++)
	{
		glm::vec2& position = positions[i];
		unsigned int particleID = particleIDs[i];

		unsigned int key = computeKey(position);

		grid[key]--;
		sortedParticleIDs[grid[key]] = particleID;
	}

}

// hash: compute hash key from the position of a particle
int HashTable::tableHash(glm::vec2& cellPos)
{
	int x = cellPos.x;
	int y = cellPos.y;

	int n = numberOfCells;

	return abs((x * P1) + ( y * P2)) % n;

}

// getCellPosition: convert world position to cell position
glm::vec2 HashTable::getCellPosition (glm::vec2& position)
{
	// discretize particle position from float to int
	int x = floor(position.x / settings->radiusOfInfluence);
	int y = floor(position.y / settings->radiusOfInfluence);
	return glm::vec2(x, y);
}

int HashTable::computeKey(glm::vec2& position)
{
	glm::vec2 cellPos = getCellPosition(position);
	unsigned int key = tableHash(cellPos);

	return key;
}