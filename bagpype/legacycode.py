    def hydrophobic_selection_sparse(self, weights, rows, columns):

        D = scipy.sparse.csr_matrix(
            (weights, (rows, columns)),
            shape=(len(self.protein.atoms), len(self.protein.atoms)),
            dtype=np.float,
        )
        D = D + D.transpose()
        node_list = (np.where(D.getnnz(0) > 0)[0]).tolist()
        _node_dictionary = dict(zip(range(len(self.protein.atoms)), node_list))

        D = D[D.getnnz(1) > 0][:, D.getnnz(0) > 0]

        N = D.shape[0]

        _Emst, LLink = self.RMST_prim_algorithm_sparse(D)
        LLink = scipy.sparse.csr_matrix(LLink)

        LLink.setdiag(-10.0)

        Dtemp = D.copy()
        # Dtemp.setdiag(D.max())

        # Dtemp[Dtemp==0] = np.amax(D)
        mD = Dtemp[np.arange(N), Dtemp.argmin(axis=0)]
        # print(mD.T.shape,  scipy.sparse.csr_matrix(np.ones([1, N])).shape  )
        mD = scipy.sparse.csr_matrix(np.ones([N, 1])) * mD
        mD = abs(mD + mD.transpose()) / 2.0
        mD *= self.hydrophobic_RMST_gamma

        E_criterion = (LLink + mD > D).astype(int)

        E_final = E_criterion.multiply(D)
        nonzeros = scipy.sparse.triu(E_final).nonzero()
        matches = list(zip(nonzeros[0].tolist(), nonzeros[1].tolist()))

        # print("    RMST sparsification used. Accepted: " + str(len(matches)-np.count_nonzero(np.triu(Emst))) + ", rejected: " + str(len(graph.edges) - len(matches)) +
        #       "; MST size: " + str(np.count_nonzero(np.triu(Emst))))

        return matches

    def RMST_prim_algorithm_sparse(self, D):
        # Initialise matrix containing all largest links along MST
        LLink = np.zeros(D.shape)
        # LLink = scipy.sparse.lil_matrix(D.shape)
        # LLink[:] = 0.0000001

        # Number of nodes in the network
        N = D.shape[0]

        # Allocate a matrix for the edge list
        E = scipy.sparse.lil_matrix(D.shape)

        # Start with a node
        mstidx = np.array([0])
        otheridx = np.arange(1, N, 1)

        T = D[otheridx, 0]
        P = np.zeros(otheridx.size, dtype=int)

        while T.size > 0:
            i = T.argmin()
            idx = otheridx[i]

            # Start with a node
            E[idx, P[i]] = D[idx, P[i]]
            E[P[i], idx] = D[P[i], idx]

            # 1) Update the longest links
            # indexes of the nodes without the parent
            tempmstidx = mstidx[mstidx != P[i]]

            # 2) update the link to the parent
            LLink[idx, P[i]] = -1 * D[idx, P[i]]
            LLink[P[i], idx] = -1 * D[P[i], idx]

            # 3) find the maximal
            if len(tempmstidx) > 0:
                new_edge = -1 * D[idx, P[i]]
                tempLLink = LLink[P[i], tempmstidx]
                tempLLink[tempLLink > new_edge] = new_edge

                LLink[idx, tempmstidx] = tempLLink
                LLink[tempmstidx, idx] = tempLLink

            # As a node is added clear his entries
            delete_row_csr(T, i)
            P = np.delete(P, i)

            # Add the node to the list
            mstidx = np.append(mstidx, idx)

            # Remove the node from the list of the free nodes
            otheridx = np.delete(otheridx, i)

            # update the distance matrix
            Ttemp = D[otheridx, idx]

            if T.size > 0:
                P[(Ttemp < T).nonzero()[0]] = idx
                T = scipy.sparse.lil_matrix.minimum(T, Ttemp)

        return E, -LLink
