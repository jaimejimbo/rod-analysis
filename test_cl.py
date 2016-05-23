import methods as m

a = [1,2,3,4,5]
b = [1,1,1,1,1]

print m.sum_arrays_with_cl(a,b)
print m.array_x_scalar_cl(a, 5)
print m.array_x_scalar_cl(a, 1.0/(sum(a)*1.0/len(a)))
