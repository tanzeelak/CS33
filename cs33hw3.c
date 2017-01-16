transpose(int *dst, int *src, int dim)
{
	int i, j, x, y, xDimen;
	for (i = 0; i < dim; i+=64)
		for (j = 0; j < dim; j+=64)
			for (x = i; x < min(i+64, dim); x++)
			{
				xDimen = x*dim;
				for (y = j; y < min(j+64, dim); y++)
					dst[x+y*dim] = src[y+xDimen];
			}
}
