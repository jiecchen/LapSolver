public class Grid2 {

    public int[] src;
    public int[] dst;
    public int[] weight;

    public Grid2(int width, int height, int clmweight) {

        //number of edges, vertices, non-bottom row
        int ne = (2 * width * height) - width - height;
        int nv = width * height;
        int nonbot = width * (height - 1);

        //i -> from, j->to, v->weight
        int src[] = new int[ne];
        int dst[] = new int[ne];
        int weight[] = new int[ne];

        //current edge
        int e = 0;

        for(int i = 1; i <= nv; i++) {
            //horizontal edges
            if ((i % width) != 0) {
                src[e] = i + 1;
                dst[e] = i;
                weight[e] = 1;
                e++;
            }

            //vertical edges
            if (i <= nonbot) {
                src[e] = i + width;
                dst[e] = i;
                weight[e] = clmweight;
                e++;
            }
        }
    }
}
