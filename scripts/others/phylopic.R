 
library(rphylopic)

setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/phyloPic")

ids = list(Discoba_Diplonema_papillatum="e0d5529a-2984-4ed1-b22b-a284d66f442d",
		   Haptista_Coccolithus_pelagicus="2a611587-6f4d-437e-a79e-4eab8fb545d3",
			Archaeplastida_Rhodophyta_Chondrus_crispus="e9df48fe-68ea-419e-b9df-441e0b208335",
			Archaeplastida_Embryophyta_Juniperus_virginiana="fbe23759-930e-43c6-a138-9d052d0a8bda",
			Archaeplastida_Embryophyta_Quercus_robur="43afe2df-ab6c-47e6-a105-c0a82b8af1c5",
			Archaeplastida_Chlorophyta_Chlorodendrophyceae="86076903-2f07-46fa-bbf8-88ba39162e63",
			Archaeplastida_Chlorophyta="e58d9dd8-6829-46f1-b49b-ed74bd780543",
			Rhizaria_Cercozoa="3580adc5-435d-4fab-a8fd-87c2f9f85f21",
			Rhizaria_Cercozoa_Paulinella_chromatophora="748f8ef9-6769-4167-930f-f8866d6a5d6f",
			Rhizaria_Endomyxa_Gromia_oviformis="ef6c6caf-12a4-4293-8fe8-759cb6ab930d",
			Rhizaria_Foraminifera_Globigerinoides_ruber="57ac4823-613c-4525-a935-41fedab41ac0",
			Rhizaria_Radiolaria_Rhizosphaera_trigonacantha="c07464fc-9719-4b16-a82b-85154f659610",
			Stramenopila_Sagenista_labyrinthulomycetes="7b10bc85-26ff-4a5d-b466-cc1269cacbf3",
			Stramenopila_Gyrista_Dictyochophyceae="a9b2c5eb-562e-4460-80d4-8898b077dc04",
			Stramenopila_Gyrista_Diatomea_Triceratium_robertsianum="4924b6bd-cfb8-4d60-a32a-442d02afbe85",
			Stramenopila_Gyrista_Saccorhiza_polyschides="da150717-f583-4c09-91e1-ffeba0744242",
			Alveolata_Ciliophora_Paramecium_tetraurelia="595359bd-6638-4900-8536-82b2ed511a25",
			Alveolata_Apicomplexa_Plasmodium_falciparum="e6014244-4dd5-4785-bf2e-c67dc4d05ca8",
			Alveolata_Dinoflagellata_Ceratium_hirundinella="573efc92-2521-453c-9fa3-0b356d78b910",
			Breviata_anathema="e789d31d-cadf-42be-9508-58be6e2a5fff",
			Amoebozoa_Amoeba_proteus="4227f7b5-b1e5-4af6-99cc-af68795f5855",
			Nucletmycea_Nucleariida_Nuclearia_moebiusi="ae4dd09d-3ffe-49fb-9054-09fdedf764fb",
			Nucletmycea_Chytridiomycota_Siphonaria_petersenii="abb669df-1335-4072-a69f-95824bf98f3c",
			Nucletmycea_Glomeromycota_Glomus_diaphanum="31aa0781-ca28-4362-8955-520f1c45f232",
			Nucletmycea_Basidiomycota_Boletus_edulis="e5d32221-7ea9-46ed-8e0a-d9dbddab0b4a",
			Nucletmycea_Ascomycota_Aspergillus_nidulans="7ebbf05d-2084-4204-ad4c-2c0d6cbcdde1",
			Holozoa_Choanoflagellata_Acanthoecidae="93f26e0e-5b4a-41a8-ae85-0f6b83d37030",
			Holozoa_Cnidaria_Medusozoa="839b9df7-c97f-444a-a611-95609e950960",
			Holozoa_Craniata_Salmo_trutta="23a7d09d-4a4d-4ad5-ad07-49a6b59a7fba",
			Holozoa_Mollusca_Myosotella_myosotis="c5835123-e2d3-4c20-9e7a-f7b6528bbf8e",
			Holozoa_Nematoda_Oscheius_dolichura="ac5299df-b34c-471a-8897-a7c10733b055",
			Holozoa_Arthropoda_Crustacea_Penaeus="aed2513d-2386-4218-b913-384838c0107b")

pics = list()
for(i in 1:length(ids)){
	cat("\r  Downloading (", i, "/", length(ids), ") ", names(ids)[i], "                    ", sep="", end="")
	# uuid = get_uuid(name = ids[[i]], n = 1)
	img = get_phylopic(uuid = ids[[i]])
	pics[i] = img
	save_phylopic(img = img, path=paste0(names(ids)[i], ".svg"))
}; rm(i, img); cat("\n")

# plot(x = 1, y = 1, type = "n")
# add_phylopic_base(img = pics[[12]], x = 1.25, y = 1.25, ysize = 0.25)
