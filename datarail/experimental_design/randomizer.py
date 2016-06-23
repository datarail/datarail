def randomizer(drug_panel, all_treatments, cntrl_pos):
    
    nwells = cntrl_pos.size   
    cntrl_idx = np.nonzero(cntrl_pos.reshape(1, nwells))[1]
    idx = np.random.choice(range(nwells),
                                 size=nwells,
                                 replace=False)
    idx[cntrl_idx] = nwells
    # idx_sort = [idx.index[x] for x in sorted(idx)]
    idx_sort = np.argsort(idx)

    panel = drug_panel.reshape(1, nwells)
    panel[0, idx_sort[:all_treatments.shape[1]]] = all_treatments[0, :]
    panel = panel.reshape(drug_panel.shape)

    return panel

          
    
                           
    
