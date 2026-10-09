ATOMS = ("rtr", "rtx", "rty", "rts", "xtr", "xtx", "xty", "xts",
         "ytr", "ytx", "yty", "yts", "str", "stx", "sty", "sts")
GROUPS = ("CN", "MCI", "AD")
PAIRS = (("CN", "MCI"), ("CN", "AD"), ("MCI", "AD"))
GROUP_COLORS = {"CN": "#f09800", "MCI": "#f00030", "AD": "#5c95ff",
                "28-30": "#f09800", "29-30": "#f09800", "30": "#f09800",
                "20-27": "#f00030", "20-28": "#f00030", "20-29": "#f00030", "0-19": "#5c95ff"}
NETWORKS = {1: "Visual", 2: "Somatomotor", 3: "Dorsal Attention", 4: "Salience / VA",
            5: "Limbic", 6: "Control", 7: "Default", 8: "Subcortical"}
NETWORK_COLORS = {1: "#f5544e", 2: "#f5c84e", 3: "#d1f54e", 4: "#4ef5cb",
                  5: "#4eadf5", 6: "#964ef5", 7: "#f54eea", 8: "gray"}
TAXA_SIX = {"upward_causation": (7, 11, 3), "downward_causation": (13, 14, 12),
            "erasure": (1, 2), "copy": (4, 8), "transfer": (6, 9), "storage": (0, 5, 10, 15)}
TAXA_FOUR = {"causation": (7, 11, 3, 13, 14, 12), "copy_erasure": (4, 8, 1, 2),
             "transfer": (6, 9), "storage": (0, 5, 10, 15)}
TAXA_PRINCIPAL = {name: TAXA_SIX[name] for name in
                  ("upward_causation", "downward_causation", "transfer", "storage")}
MOCA_GROUPS = ("29-30", "20-28", "0-19")
MMSE_GROUPS = ("28-30", "20-27", "0-19")
