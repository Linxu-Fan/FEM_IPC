#include "objMesh.h"



void objMeshFormat::readObjFile(std::string fileName)
{
	clear();

	// read vertices and faces
	std::ifstream in;
	in.open(fileName);
	std::string line;
	while (getline(in, line))
	{
		if (line.size() > 0)
		{
			std::vector<std::string> vecCoor = split(line, " ");
			if (vecCoor[0] == "v")
			{
				Eigen::Vector3d vt = { std::stod(vecCoor[1]) , std::stod(vecCoor[2]) , std::stod(vecCoor[3])};
				vertices.push_back(vt);
			}

			if (vecCoor[0] == "f")
			{
				Eigen::Vector3i ele = { std::stoi(vecCoor[1]) - 1 ,  std::stoi(vecCoor[2]) - 1 ,  std::stoi(vecCoor[3]) - 1 };
				faces.push_back(ele);
			}

		}
	}
	in.close();



	findVertFaces_Edges();

}

void objMeshFormat::findVertFaces_Edges()
{
	std::set<std::string> edgesStored;
	for (int i = 0; i < faces.size(); i++)
	{
		int v1 = faces[i][0], v2 = faces[i][1], v3 = faces[i][2];
		std::string v1v2 = std::to_string(std::min(v1, v2)) + std::to_string(std::max(v1, v2));
		std::string v1v3 = std::to_string(std::min(v1, v3)) + std::to_string(std::max(v1, v3));
		std::string v2v3 = std::to_string(std::min(v2, v3)) + std::to_string(std::max(v2, v3));
		if (edgesStored.find(v1v2) == edgesStored.end())
		{
			edgesStored.insert(v1v2);
			Eigen::Vector2i eg = { std::min(v1, v2) ,  std::max(v1, v2) };
			edges.push_back(eg);
		}

		if (edgesStored.find(v1v3) == edgesStored.end())
		{
			edgesStored.insert(v1v3);
			Eigen::Vector2i eg = { std::min(v1, v3) ,  std::max(v1, v3) };
			edges.push_back(eg);
		}

		if (edgesStored.find(v2v3) == edgesStored.end())
		{
			edgesStored.insert(v2v3);
			Eigen::Vector2i eg = { std::min(v2, v3) ,  std::max(v2, v3) };
			edges.push_back(eg);
		}
	}


	// find vertFaces
	vertFaces.resize(vertices.size());
	for (int fi = 0; fi < faces.size(); fi++)
	{
		Eigen::Vector3i fc = faces[fi];
		vertFaces[fc[0]].push_back(fi);
		vertFaces[fc[1]].push_back(fi);
		vertFaces[fc[2]].push_back(fi);
	}


}

void depthFirstSearch(int v, const std::vector<std::set<int>>& adjList, std::vector<bool>& visited, std::vector<int>& component) 
{
	visited[v] = true;
	component.push_back(v);
	for (int neighbor : adjList[v]) 
	{
		if (!visited[neighbor]) 
		{
			depthFirstSearch(neighbor, adjList, visited, component);
		}
	}
}


void objMeshFormat::sepConnectedComponents()
{

	// find connected components
	std::vector<std::vector<int>> components;
	std::vector<std::set<int>> adjList(vertices.size());
	for (const auto& f : faces) 
	{
		adjList[f[0]].insert(f[1]);
		adjList[f[0]].insert(f[2]);
		adjList[f[1]].insert(f[0]);
		adjList[f[1]].insert(f[2]);
		adjList[f[2]].insert(f[0]);
		adjList[f[2]].insert(f[1]);
	}
	std::vector<bool> visited(vertices.size(), false);
	for (size_t i = 0; i < vertices.size(); ++i) 
	{
		if (!visited[i]) 
		{
			std::vector<int> component;
			depthFirstSearch(i, adjList, visited, component);
			components.push_back(component);
		}
	}

	// remove unnecessary vertices
	for (size_t i = 0; i < components.size(); ++i) 
	{
		objMeshFormat cop;
		std::unordered_map<int, int> vertexMap;

		for (int idx : components[i]) {
			vertexMap[idx] = cop.vertices.size();
			cop.vertices.push_back(vertices[idx]);
		}

		for (const auto& f : faces) 
		{
			if (vertexMap.count(f[0]) && vertexMap.count(f[1]) && vertexMap.count(f[2]))
			{
				cop.faces.push_back({ vertexMap[f[0]], vertexMap[f[1]], vertexMap[f[2]] });
			}
		}
		
		cop.findVertFaces_Edges();
		componentsSep.push_back(cop);
	}

}



void objMeshFormat::clear()
{
	vertices.clear();
	faces.clear();
	vertFaces.clear();
	componentsSep.clear();
}


void objMeshFormat::outputFile(std::string fileName, int timestep)
{
	std::ofstream outfile9("./output/" + fileName + "_" + std::to_string(timestep) + ".obj", std::ios::trunc);
	for (int k = 0; k < vertices.size(); k++)
	{
		Eigen::Vector3d scale = vertices[k];
		outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
	}
	for (int k = 0; k < faces.size(); k++)
	{
		outfile9 << "f ";
		if (0)
		{
			for (int m = 0; m < faces[k].size(); m++)
			{
				outfile9 << faces[k][m] + 1 << " ";
			}
		}
		else
		{
			for (int m = faces[k].size() - 1; m >= 0; m--)
			{
				outfile9 << faces[k][m] + 1 << " ";
			}
		}
		outfile9 << std::endl;
	}
	outfile9.close();
}


