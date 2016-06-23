def randomizer(drug_panels, all_treatments, cntrl_pos):
    
    nwells = cntrl_pos.size   
    cntrl_idx = np.nonzero(cntrl_pos.reshape(1, nwells))[1]
    idx = np.random.choice(range(nwells),
                                 size=nwells,
                                 replace=False)
    idx[cntrl_idx] = nwells
    idx_sort = np.argsort(idx)
    panels = []
    for i, panel in enumerate(drug_panels):
        panel = drug_panel.reshape(1, nwells)
        panel[0, idx_sort[:all_treatments.shape[1]]] = all_treatments[i, :]
        panel = panel.reshape(drug_panel.shape)
        panels.append(panel)
    return panels

          
    
                           
    
