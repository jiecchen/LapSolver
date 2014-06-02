#ifndef PROFILER_HPP
#define PROFILER_HPP

#include "graph.hpp"
#include <stdio.h>
#include <time.h>
#include <vector>

class Profiler
{
public:
    Profiler(bool disp)
	: display(disp), last_log(0.0), start_time(0.0), end_time(0.0), cluster_buffer(-1, 0.0)
	{

	}
    ~Profiler() {}
    void log(const char* event)
	{
	    double cur_time = get_time();
	    if (display)
		printf("%lf seconds (%lf): %s\n", cur_time, cur_time - last_log, event);
	    last_log = cur_time;
	}
    void new_cluster(int size, double conductance)
	{
	    if (cluster_buffer.first > 0)
		clusters.push_back(cluster_buffer);
//	    printf("size: %d, conduct: %lf\n", size, conductance);
	    cluster_buffer.first = size;
//	    if (size == 1)
//		cluster_buffer.second = 0.0;
//	    else
		cluster_buffer.second = conductance;
	}
    void new_recluster(int size, double conductance)
	{
	    if (cluster_buffer.first > 0)
		rejected_clusters.push_back(cluster_buffer);
	    else
		exit(666);  // should already have a cluster if reclustering!

	    cluster_buffer.first = size;
//	    if (size == 1)
//		cluster_buffer.second = 0.0;
//	    else
		cluster_buffer.second = conductance;
	}

    // only call after end_alg!
    void cluster_stats()
	{
	    std::vector<double> means;
	    means = cluster_mean(&clusters);
	    printf("Ended with %d clusters (mean size = %lf, mean conductance = %lf ignoring %d singles)\n", 
		   int(clusters.size()), means[0], means[1], int(means[2]));
	    means = cluster_mean(&rejected_clusters);
	    printf("Rejected %d clusters (mean size = %lf, mean conductance = %lf ignoring %d singles)\n", 
		   int(rejected_clusters.size()), means[0], means[1], int(means[2]));
	}
    void write_clust_info(FILE* fp)
	{
	    int n = clusters.size();
	    fwrite(&n, sizeof(int), 1, fp);
	    //printf("n=%d\n", n);
	    for (unsigned int x = 0; x < clusters.size(); x++)
	    {
		fwrite(&clusters[x].first, sizeof(int), 1, fp);
	    }
	    for (unsigned int x = 0; x < clusters.size(); x++)
	    {
		fwrite(&clusters[x].second, sizeof(double), 1, fp);
//		printf("%lf\n", clusters[x].second);
	    }
	    n = rejected_clusters.size();
	    fwrite(&n, sizeof(int), 1, fp);
	    for (unsigned int x = 0; x < rejected_clusters.size(); x++)
	    {
		fwrite(&rejected_clusters[x].first, sizeof(int), 1, fp);
	    }
	    for (unsigned int x = 0; x < rejected_clusters.size(); x++)
	    {
		fwrite(&rejected_clusters[x].second, sizeof(double), 1, fp);
	    }
	}
    void start_alg()
	{
	    start_time = get_time();
	}
    void end_alg()
	{
	    end_time = get_time();
	    // flush the cluster buffer
	    if (cluster_buffer.first > 0)
	    {
		clusters.push_back(cluster_buffer);
	    }
	}
    void log_alg()
	{
	    printf("***Algorithm itself took %lf seconds.***\n", end_time - start_time);
	}
private:
    bool display;
    double last_log;
    double start_time;
    double end_time;
    std::vector<std::pair<int, double> > clusters;
    std::vector<std::pair<int, double> > rejected_clusters;
    std::pair<int, double> cluster_buffer; // when a new cluster is reported buffer it here

    std::vector<double> cluster_mean(std::vector<std::pair<int, double> >* vec)
	{
	    double mean_clust = 0;
	    double mean_conduct = 0;
	    int num_single_clusts = 0;
	    for (unsigned int x = 0; x < vec->size(); x++)
	    {
		mean_clust += (*vec)[x].first;
		if ((*vec)[x].first > 1)
		    mean_conduct += (*vec)[x].second;
		else
		    num_single_clusts++;
	    }
	    mean_clust /= (double)vec->size();
	    mean_conduct /= (double(vec->size()) - double(num_single_clusts));
	    std::vector<double> out;
	    out.push_back(mean_clust);
	    out.push_back(mean_conduct);
	    out.push_back(num_single_clusts);
	    return out;
	}
    double get_time()
	{
	    return (double)clock()/(double)CLOCKS_PER_SEC;
	}
};


#endif // PROFILER_HPP
