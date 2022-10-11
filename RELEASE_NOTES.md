## 3.0.0
    - BACKWARDS INCOMPATIBILIY: the output structure for the `save_assemblies_from_fastas`
      method has changed.
    - DEPRECATION: `save_assembly_from_fasta`.
    - Added the `save_assembly_from_fasta2` function.
    - Updated `save_assemblies_from_fastas` to return the path to the filtered input file as
      well as the workspace UPA.

## 2.0.1
### Update
    - No changes, just a version bump to allow registering the most recent version of the code on
      the CI environment.

## 2.0.0
### Update
    - Added the save_assemblies_from_fastas_function
    - Removed the unused UI import apps.
    - Removed the `ftp_url` parameter from `save_assembly_from_fasta`
    - Removed the `taxon_ref` parameter from `save_assembly_from_fasta`. It is now silently
      ignored.
    - Clarified `save_assembly_from_fasta` documentation
    - Fixed a bug in `save_assembly_from_fasta` that could cause file name collisions if
      multiple `AssemblyUtil` instances are run in parallel
    - Clarified the error message text when an empty file is submitted or results from
      contig filtering, and prevented a guaranteed to fail attempt to save the file to the
      Blobstore

## 1.2.6
### Update
	- Added unit test for error message

## 1.2.5
### Update
	- Change error message for uploading empty message

## 1.2.4
### Fixed
	- Bug fix for test objects

## 1.2.3
### Fixed
	- Bug fix for export_as_fasta

## 1.2.2:
### Feature
	- adding AnnotatedMetagenomeAssembly functionality to get_fastas function

## 1.2.1:
### Feature
	- Added get_fastas function in AssemblyUtilImpl

## 1.1.0:
### Update
	- Update to Python 3