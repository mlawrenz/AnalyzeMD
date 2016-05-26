import copy
import collections

def keypaths(nested):
    for key, value in nested.iteritems():
        if isinstance(value, collections.Mapping):
            for subkey, subvalue in keypaths(value):
                yield [key] + subkey, subvalue
        else:
            yield [key], value

example=dict()
example['A']=dict()
example['B']=dict()
example['A']['a1']=1
example['A']['a2']=2
example['B']['b2']=3
example['B']['b1']=4

reverse_dict = {value: keypath for keypath, value in keypaths(example)}
print reverse_dict
chain=reverse_dict[4][0]
res=reverse_dict[4][1]




