# Partial documented liquid ledger. Values are volumes and analytical moles.
# No acid speciation, density, mixing enthalpy, or full waste prediction.
v=Me.InputVariables
o=Me.OutputVariables
o.Clear()
total=0.0
hno3=0.0
hcl=0.0
for i in range(int(v['recipe_rows'])):
    prefix='row_'+str(i)+'_'
    vol=v[prefix+'mL']
    total+=vol
    hno3+=vol*v[prefix+'HNO3_M']/1000.0
    hcl+=vol*v[prefix+'HCl_M']/1000.0
o['documented_fresh_volume_mL']=total
o['documented_internal_transfer_mL']=v['internal_transfer_mL']
o['documented_pumped_volume_mL']=total+v['internal_transfer_mL']
o['documented_pumping_time_min']=(total+v['internal_transfer_mL'])/v['flow_mL_min']
o['unallocated_time_min']=v['reported_process_min']-o['documented_pumping_time_min']
o['analytical_HNO3_mol']=hno3
o['analytical_HCl_mol']=hcl
o['unallocated_fraction']=o['unallocated_time_min']/v['reported_process_min']
o['partial_nonproduct_collection_mL']=total-v['final_collection_mL']
for _name in list(o.Keys):
    globals()[str(_name)] = o[_name]
