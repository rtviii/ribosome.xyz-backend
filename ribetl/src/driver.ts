import path from 'path'
import fs from 'fs'
import yargs from 'yargs';
import shell from 'shelljs'
import {processPDBRecord} from './requestGqlProfile'
import { RibosomeStructure } from './RibosomeTypes';
import * as dotenv from "dotenv";

dotenv.config({ path: '/home/rxz/dev/riboxyzbackend/rxz_backend/.env' });

export const writeupdateStruct = (r:RibosomeStructure) =>{

    var   rcsb_id         = r.rcsb_id
    var   target_filename = path.join(process.env.STATIC_ROOT as string, rcsb_id.toUpperCase(),rcsb_id.toUpperCase() +'.json')

    if (!fs.existsSync(path.dirname(target_filename))){
        shell.mkdir('-p', path.dirname(target_filename))
    }

      fs.writeFileSync(target_filename, JSON.stringify(r, null, 4))
      console.log("Has written to ", path.dirname(target_filename))

}

export const updateCatalogueWStruct = async (
  new_entry:CatalogueEntry) =>{

    if (fs.existsSync(process.env.CATALOGUE as string)){

    } else{
      console.log("Catalogue file doesn't exist. Creating.");
      fs.writeFileSync(process.env.CATALOGUE as string,JSON.stringify({}))
    }
        var catalogue = JSON.parse(fs.readFileSync(process.env.CATALOGUE as string, 'utf-8'))
        var rcsb_id:string = Object.keys(new_entry)[0]
        if ( Object.keys(catalogue).includes(rcsb_id ) ){
          catalogue[rcsb_id] = new_entry
        }
        else{
        catalogue = {...catalogue, ...new_entry}
        }
        fs.writeFileSync(process.env.CATALOGUE as string, JSON.stringify(catalogue, null, 4))
    }

(async () =>{


	const args = await yargs.options(
    {'struct': {type:'string',demandOption:true, alias:"s" }}
    ).argv

  if (args.struct) {
    if (args.struct.length < 4){
      console.log("Enter a valid RCSB PDB ID to build from RCSB's gql api.");
      process.exit(2)
    }
    console.log(`Processing ${args.struct}`);
    
    try{
      var struct: RibosomeStructure = await processPDBRecord(args.struct)
      var new_entry: CatalogueEntry = {
        [struct.rcsb_id.toUpperCase()]: {
          _did_rename_chains: false,
          _did_split_chains: false,
          _ligands_parsed: false,
          status: true,
          error: null
        }
      };
      await updateCatalogueWStruct(new_entry)
      writeupdateStruct(struct)
    }
    catch (e: any) {
      var new_entry: CatalogueEntry = {
        [args.struct.toUpperCase()]: {
          _did_rename_chains: false,
          _did_split_chains: false,
          _ligands_parsed: false,
          status: true,
          error: e
        }
      };
	  console.log("Failed: \n\n");
	  console.log("+++++++++++++++++")
    console.log(e);
    console.log("+++++++++++++++++")

      
    await updateCatalogueWStruct(new_entry)

    }
  }
})()




type CatalogueEntry = {
  [rcsb_id: string]: {
    status            : boolean,
    _ligands_parsed   : boolean,
    _did_rename_chains: boolean
    _did_split_chains : boolean,
    error: string | null
  }
}

