#!/usr/bin/env python
# dx-instance-type.py
# outputs the available instances of a region that a particular job is located
import dxpy
import argparse
from subprocess import check_call
from dxpy.utils.completer import InstanceTypesCompleter


def parse_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Get a list of known instance types')

    ap.add_argument('-m', '--min-memory',
                    help='Minimum memory in GB.',
                    type=int,
                    required=False)

    ap.add_argument('-c', '--min-cores',
                    help='Minimum number of cores.',
                    type=int,
                    required=False)

    ap.add_argument('-s', '--min-storage',
                    help='Minimum storage in GB.',
                    type=int,
                    required=False)

    ap.add_argument('-p', '--project',
                    help='Restrict to instances valid for current project.',
                    action='store_true',
                    required=False)

    args = ap.parse_args()

    return args


def make_instance_type_to_resources():
    """ Returns a dictionary of keys-instance type and values-resource dictionary """
    instance_type_to_resources = {it.Name: {'memory_gb': it.Memory_GB, 'storage_gb': it.Storage_GB, 'cores': it.CPU_Cores}
        for it in InstanceTypesCompleter.instance_types.values()}

    return instance_type_to_resources


def filter_instances(instance_type_to_resources, region):
    """ filters the instance dictionary to the available instances on a given cloud service """

    if 'azure' in region:
        filtered_dict = {k: v for k, v in instance_type_to_resources.iteritems() if 'azure' in k}
    elif region == 'aws:us-east-1':
        filtered_dict = {k: v for k, v in instance_type_to_resources.iteritems() if 'azure' not in k}
    elif 'aws' in region:
        # AWS regions outside of aws:us-east-1 do not offer hdd-based instances.
        filtered_dict = {k: v for k, v in instance_type_to_resources.iteritems() if 'azure' not in k and 'hdd' not in k}
    else:
        # @TODO: Should we print a warning if we don't recognize this region?
        filtered_dict = instance_type_to_resources,

    return filtered_dict



def output_region_compatible_instances(instance_type_to_resources):
    def instance_compare(a):
        mem, storage, cores = a.split('_')
        return (mem, storage, int(cores[1:]))

    sorted_keys = sorted(instance_type_to_resources.keys(), key=instance_compare)
    #print 'Instance Type\tMemory (GB)\tStorage (GB)\tCores'
    print 'Name\tMemory_GB\tStorage_GB\tCPU_Cores'
    for inst_type in sorted_keys:
        resource_dict = instance_type_to_resources[inst_type]
        print '{}\t{}\t{}\t{}'.format(inst_type, resource_dict['memory_gb'], resource_dict['storage_gb'], resource_dict['cores'])


def main(min_memory, min_cores, min_storage, restrict_to_current_project):
    instance_type_to_resources = make_instance_type_to_resources()

    # If we are in a job, then filter the instances to only those available in
    # the region for the given job.
    if restrict_to_current_project:
        instance_type_to_resources = filter_instances(instance_type_to_resources,
            dxpy.api.job_describe(dxpy.PROJECT_CONTEXT_ID)['region'])

    if min_memory is not None:
        instance_type_to_resources = {k: v for k, v in instance_type_to_resources.iteritems() if v['memory_gb'] >= min_memory}
    if min_storage is not None:
        instance_type_to_resources = {k: v for k, v in instance_type_to_resources.iteritems() if v['storage_gb'] >= min_storage}
    if min_cores is not None:
        instance_type_to_resources = {k: v for k, v in instance_type_to_resources.iteritems() if v['cores'] >= min_cores}

    output_region_compatible_instances(instance_type_to_resources)

if __name__=="__main__":
    args = parse_args()
    main(args.min_memory, args.min_cores, args.min_storage, args.project)
