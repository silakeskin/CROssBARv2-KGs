from __future__ import annotations

import collections
import csv
import re
import asyncio
import warnings

from concurrent.futures.thread import ThreadPoolExecutor

from abc import ABC, abstractmethod
from collections.abc import Iterable

import pypath.resources.urls as urls
from pypath.share import curl

import asyncio
import aiohttp
import warnings
import random

from dataclasses import dataclass
from typing import Optional
import json




DEFAULT_TIMEOUT = aiohttp.ClientTimeout(total=30)
DDI_BATCH_SIZE = 500

_url = urls.urls['kegg_api']['url']

class AsyncRateLimiter:
    """
    Ensures that requests are started no faster than the configured rate.
    """

    def __init__(self, requests_per_second: float = 2.0):
        if requests_per_second <= 0:
            raise ValueError(
                "requests_per_second must be greater than zero."
            )

        self.minimum_interval = 1.0 / requests_per_second
        self._lock = asyncio.Lock()
        self._last_request_time = 0.0

    async def wait(self) -> None:
        async with self._lock:
            loop = asyncio.get_running_loop()
            current_time = loop.time()

            elapsed = current_time - self._last_request_time
            remaining = self.minimum_interval - elapsed

            if remaining > 0:
                await asyncio.sleep(remaining)

            self._last_request_time = loop.time()

def gene_to_pathway(org):


    return _kegg_from_source_to_target('gene', 'pathway', org)


def pathway_to_gene(org):

    return _kegg_from_source_to_target('pathway', 'gene', org)


def gene_to_drug(org):

    return _kegg_from_source_to_target('gene', 'drug', org)


def drug_to_gene(org):

    return _kegg_from_source_to_target('drug', 'gene', org)


def gene_to_disease(org):

    return _kegg_from_source_to_target('gene', 'disease', org)


def disease_to_gene(org):

    return _kegg_from_source_to_target('disease', 'gene', org)


def pathway_to_drug():

    return _kegg_from_source_to_target('pathway', 'drug')


def drug_to_pathway():

    return _kegg_from_source_to_target('drug', 'pathway')


def pathway_to_disease():
    
    return _kegg_from_source_to_target('pathway', 'disease')


def disease_to_pathway():

    return _kegg_from_source_to_target('disease', 'pathway')


def disease_to_drug():

    return _kegg_from_source_to_target('disease', 'drug')


def drug_to_disease():
    
    return _kegg_from_source_to_target('drug', 'disease')


def drug_to_drug(
    drugs: list | tuple =None,
    join: bool=True, 
    asynchronous: bool=False
) -> tuple:
    """
    Downloads drug-drug interaction data from KEGG database.

    Arguments:
    @drugs: Drug IDs as a list or a tuple.
    @join: if it's True, returns individual interactions of queried list.
            Else, joins them together and returns mutual interactions.
    @asynchronous: This function yet to be implemented.
    """

    DrugToDrugInteraction = collections.namedtuple(
        'DrugToDrugInteraction', 
        (
            'type',
            'name',
            'interactions',
        ),
    )

    drug = _Drug()
    compound = _Compound()

    
    if drugs != None:

        entries = _kegg_ddi(drugs, join=join, asynchronous=asynchronous)

    else:
        drugIds = drug.get_data().keys()
        entries = _kegg_ddi(drugIds, join=False, asynchronous=asynchronous)

    interactions = dict()

    for entry in entries:

        for i in range(2):

            try:
                entry_type, entry_id = entry[i].split(':')
            except ValueError:
                entry_type = entry[i][0]
                entry_id = entry[i]

            if entry_type == 'dr' or entry_type == 'D':

                entry_type = 'drug'
                entry_db = drug

            elif entry_type == 'cpd' or entry_type == 'C':

                entry_type = 'compound'
                entry_db = compound
            
            else:

                print(f'Unknown type \'{entry_type}\', exiting...')
                exit()

            try:

                entry_name = entry_db.get(entry_id)

            except:

                entry_name = None

            tmp_dict = {
                'type': entry_type,
                'id': entry_id,
                'name': entry_name
            }
            
            if i == 0:

                source = tmp_dict
            else:
                
                target = tmp_dict
        
        labels = entry[2].split(',')
        contraindication = True if 'CI' in labels else False
        precaution = True if 'P' in labels else False

        Interaction = collections.namedtuple(
                f'{target["type"].capitalize()}Interaction',
                (
                    'type',
                    'id',
                    'name',
                    'contraindication',
                    'precaution',
                )
        )

        interaction = Interaction(
                    target['type'],
                    target['id'],
                    target['name'],
                    contraindication,
                    precaution,
        )

        diseaseId = source['id']

        try:
            interactions[diseaseId]['interactions'].append(interaction)

        except KeyError:
            interactions[diseaseId] = dict()
            interactions[diseaseId]['type'] = source['type']
            interactions[diseaseId]['name'] = source['name']
            interactions[diseaseId]['interactions'] = [interaction]

    for key, value in interactions.items():

        interactions[key] = DrugToDrugInteraction(
            value['type'],
            value['name'],
            tuple(value['interactions']),
        )
    
    return interactions

def get_diseases(diseases):

    result = _kegg_get(diseases)

    DiseaseEntry = collections.namedtuple(
        'DiseaseEntry',
        (
            'db_links',
            'references',
        )
    )

    entries = list()

    db_links = dict()
    references = list()

    db_links_regex = r'(?:DBLINKS)?\s*([^:\s]+)\s*:\s*(.+)'
    db_links_matcher = re.compile(db_links_regex)

    references_regex = r'REFERENCE\s*([^\s]+)'
    references_matcher = re.compile(references_regex)

    state = None

    for line in result:

        line = line.strip(' ')

        if line.startswith('///'):
            entries.append(
                DiseaseEntry(
                    db_links,
                    references
                )
            )

            db_links = dict()
            references = list()
            state = None
        
        elif line.startswith('DBLINKS'):
            state = 'DBLINKS'
            key, value = db_links_matcher.findall(line)[0]
            if ' ' in value:
                value = value.split(' ')
            db_links[key] = value
        
        elif line.startswith('REFERENCE'):
            state = None
            if references_matcher.findall(line):
                reference = references_matcher.findall(line)[0]
                references.append(reference)
        
        else:
            if state == 'DBLINKS':
                key, value = db_links_matcher.findall(line)[0]
                if ' ' in value:
                    value = value.split(' ')
                db_links[key] = value
            else:
                continue
    
    return entries


def kegg_gene_id_to_ncbi_gene_id(org):

    return _kegg_conv(org, 'ncbi-geneid', target_split=True)


def ncbi_gene_id_to_kegg_gene_id(org):

    return _kegg_conv('ncbi-geneid', org, source_split=True)


def kegg_gene_id_to_uniprot_id(org):

    return _kegg_conv(org, 'uniprot', target_split=True)


def uniprot_id_to_kegg_gene_id(org):

    return _kegg_conv('uniprot', org, source_split=True)


def kegg_drug_id_to_chebi_id():

    return _kegg_conv('drug', 'chebi', source_split=True, target_split=True)


def chebi_id_to_kegg_drug_id():

    return _kegg_conv('chebi', 'drug', source_split=True, target_split=True)


def _kegg_general(operation, *arguments, split=True):

    url = _url % operation

    for argument in arguments:

        url += f'/{argument}'

    c = curl.Curl(url, silent = True, large = False)

    try:
        return [line.split('\t') if split else line for line in c.result.split('\n') if line]
    except AttributeError:
        return []



async def _fetch_with_retry(
    session: aiohttp.ClientSession,
    semaphore: asyncio.Semaphore,
    rate_limiter: AsyncRateLimiter,
    operation: str,
    identifier: str,
    retries: int = 3,
):


    url = _url % operation
    url += f"/{identifier}"

    last_exception = None

    for attempt in range(1, retries + 1):

        try:

            await rate_limiter.wait()

            async with semaphore:

                async with session.get(url) as response:

                    status = response.status

                    response.raise_for_status()

                    text = await response.text()

                    return RequestResult(
                        identifier=identifier,
                        endpoint=operation,
                        success=True,
                        retries=attempt,
                        status_code=status,
                        error=None,
                        data=[
                            line.split("\t")
                            for line in text.splitlines()
                            if line
                        ],
                    )

        except aiohttp.ClientResponseError as e:

            last_exception = e

            # 404 -> DDI kaydı yok, retry yapma
            if e.status == 404:

                return RequestResult(
                    identifier=identifier,
                    endpoint=operation,
                    success=False,
                    retries=attempt,
                    status_code=e.status,
                    error=str(e),
                    data=[],
                )

            # 403 -> Rate limit olabilir, daha uzun bekle
            if e.status == 403:

                if attempt < retries:

                    retry_after = e.headers.get("Retry-After")

                    if retry_after is not None:

                        wait = float(retry_after)

                    else:

                        wait = (5 * attempt) + random.uniform(0, 2)

                    print(
                        f"[403] {identifier} "
                        f"(attempt {attempt}/{retries}) "
                        f"retrying in {wait:.1f}s"
                    )

                    await asyncio.sleep(wait)

                    continue

            # Diğer HTTP hataları (500, 502, vb.)
            if attempt < retries:

                wait = (2 ** (attempt - 1)) + random.uniform(0, 1)

                print(
                    f"[HTTP {e.status}] {identifier} "
                    f"(attempt {attempt}/{retries}) "
                    f"retrying in {wait:.1f}s"
                )

                await asyncio.sleep(wait)

        except (
            aiohttp.ClientError,
            asyncio.TimeoutError,
        ) as e:

            last_exception = e

            if attempt < retries:
                await asyncio.sleep((2 ** (attempt - 1)) + random.random())

    status = (
        getattr(last_exception, "status", None)
        if last_exception
        else None
    )

    return RequestResult(
        identifier=identifier,
        endpoint=operation,
        success=False,
        retries=retries,
        status_code=status,
        error=repr(last_exception),
        data=[],
    )

async def _kegg_general_async(
    session: aiohttp.ClientSession,
    semaphore: asyncio.Semaphore,
    rate_limiter: AsyncRateLimiter,
    operation: str,
    *arguments,
):


    if len(arguments) != 1:

        raise ValueError(
            "_kegg_general_async currently supports a single identifier."
        )

    return await _fetch_with_retry(
        session=session,
        semaphore=semaphore,
        rate_limiter=rate_limiter,
        operation=operation,
        identifier=arguments[0],
    )

def _kegg_list(database, option=None, org=None):

    if database == 'brite' and option != None:

        return _kegg_general('list', database, option)

    if database == 'pathway' and org != None:

        return _kegg_general('list', database, org)
    
    return _kegg_general('list', database)

def _kegg_get(db_entries: list | tuple | str):

    if isinstance(db_entries, list) or isinstance(db_entries, tuple):
        db_entries = '+'.join(db_entries)
    elif isinstance(db_entries, str):
        pass
    else:
        print(f'Unrecognized db_entries type: {type(db_entries)}')
        print('Exiting...')
        exit()
    
    return _kegg_general('get', db_entries, split=False)

def _kegg_conv(source_db, target_db, source_split=False, target_split=False):

    result = _kegg_general('conv', target_db, source_db)
    conversion_table = dict()
    keys = set()

    for index, entry in enumerate(result):

        source = entry[0]
        target = entry[1]

        if source_split:
            source = source.split(':')[1]

        if target_split:
            target = target.split(':')[1]

        if source in keys:
            try:
                conversion_table[source].append(target)
            except AttributeError:
                conversion_table[source] = [conversion_table[source]]
                conversion_table[source].append(target)
        else:
            conversion_table[source] = target

        keys.add(source)

    return conversion_table


def _kegg_link(source_db, target_db):

    return _kegg_general('link', target_db, source_db)


def _kegg_ddi(drugIds, join=True, asynchronous=False):

    if join and not isinstance(drugIds, str):
        drugIds = ['+'.join(drugIds)]

    if asynchronous:
        return asyncio.run(_kegg_ddi_async(drugIds))

    return _kegg_ddi_sync(drugIds)


def _kegg_ddi_sync(drugIds):

    result = list()

    if isinstance(drugIds, Iterable):

        for drugId in drugIds:

            result.extend(_kegg_general('ddi', drugId))

        return result


async def _kegg_ddi_async(
    drugIds,
    batch_size=25,
    concurrency=2,
    requests_per_second=2.0,
):
    result = []

    stats = {
        "success_first_try": 0,
        "success_after_retry": 0,
        "not_found": 0,
        "failed": 0,
    }

    status_counter = {}
    failures = []

    timeout = aiohttp.ClientTimeout(total=30)

    semaphore = asyncio.Semaphore(concurrency)

    rate_limiter = AsyncRateLimiter(
        requests_per_second=requests_per_second
    )

    async with aiohttp.ClientSession(timeout=timeout) as session:

        for batch_start in range(
            0,
            len(drugIds),
            batch_size,
        ):
            batch = drugIds[
                batch_start:
                batch_start + batch_size
            ]

            tasks = [
                _kegg_general_async(
                    session,
                    semaphore,
                    rate_limiter,
                    "ddi",
                    drug_id,
                )
                for drug_id in batch
            ]

            batch_results = await asyncio.gather(
                *tasks,
                return_exceptions=True,
            )

            for request in batch_results:

                if isinstance(request, Exception):
                    warnings.warn(
                        f"DDI request failed: {request}"
                    )

                    stats["failed"] += 1

                    status_counter["EXCEPTION"] = (
                        status_counter.get("EXCEPTION", 0) + 1
                    )

                    failures.append({
                        "drug_id": None,
                        "status_code": None,
                        "retries": None,
                        "error": repr(request),
                    })

                    continue

                print(
                    f"{request.identifier:8} | "
                    f"status={request.status_code} | "
                    f"success={request.success} | "
                    f"retries={request.retries}"
                )

                status_label = (
                    str(request.status_code)
                    if request.status_code is not None
                    else "NO_STATUS"
                )

                status_counter[status_label] = (
                    status_counter.get(status_label, 0) + 1
                )

                if request.success:
                    if request.retries == 1:
                        stats["success_first_try"] += 1
                    else:
                        stats["success_after_retry"] += 1

                    if request.data:
                        result.extend(request.data)

                    continue


                if request.status_code == 404:
                    stats["not_found"] += 1
                    continue

                stats["failed"] += 1

                failures.append({
                    "drug_id": request.identifier,
                    "status_code": request.status_code,
                    "retries": request.retries,
                    "error": request.error,
                })

            processed = min(
                batch_start + len(batch),
                len(drugIds),
            )

            print(
                f"Processed {processed}/{len(drugIds)} drugs "
                f"({processed / len(drugIds) * 100:.1f}%)"
            )

    with open(
        "ddi_failures.json",
        "w",
        encoding="utf-8",
    ) as file:
        json.dump(
            failures,
            file,
            indent=2,
            ensure_ascii=False,
        )

    print("\n========== DDI DOWNLOAD SUMMARY ==========")

    print(
        f"Succeeded on first try : "
        f"{stats['success_first_try']}"
    )

    print(
        f"Succeeded after retry  : "
        f"{stats['success_after_retry']}"
    )

    print(
        f"No DDI record (404)     : "
        f"{stats['not_found']}"
    )

    print(
        f"Failed permanently     : "
        f"{stats['failed']}"
    )

    print(
        "Failure log            : "
        "ddi_failures.json"
    )

    print("\nStatus code summary")

    for status, count in sorted(status_counter.items()):
        print(f"{status}: {count}")

    print(
        f"\nReturned interaction rows: "
        f"{len(result)}"
    )

    return result

def _kegg_from_source_to_target(source_db, target_db, org=None) -> tuple:

    db_name_list = [source_db, target_db]
    db_list = list()

    for db in db_name_list:

        if db == 'gene' and org != None:

            db_list.append(_Gene(org))

            kegg_to_ncbi = _KeggToNcbi(org)
            kegg_to_uniprot = _KeggToUniprot(org)

        elif db == 'pathway':

            db_list.append(_Pathway())

        elif db == 'disease':

            db_list.append(_Disease())

        elif db == 'drug':

            db_list.append(_Drug())

            kegg_to_chebi = _KeggToChebi()

        else:
            print('Problem in function call. Check arguments.')
            exit()
        
    if target_db == 'gene':
        TargetDbEntry = collections.namedtuple(
            f'{target_db.capitalize()}Entry',
            [
                f'{target_db}_id',
                f'{target_db}_name',
                'ncbi_gene_id',
                'uniprot_ids'
            ]
        )

    elif target_db == 'drug':
        TargetDbEntry = collections.namedtuple(
            f'{target_db.capitalize()}Entry',
            [
                f'{target_db}_id',
                f'{target_db}_name',
                'chebi_id',
            ]
        )

    else:
        TargetDbEntry = collections.namedtuple(
            f'{target_db.capitalize()}Entry',
            [
                f'{target_db}_id',
                f'{target_db}_name',
            ]
        )

    
    if source_db == 'gene':
        Interaction = collections.namedtuple(
            f'{source_db.capitalize()}To{target_db.capitalize()}Interaction',
            [
                f'{source_db}_name',
                f'{target_db.capitalize()}Entries',
                'ncbi_gene_id',
                'uniprot_ids'
            ]
        )

    elif source_db == 'drug':
        Interaction = collections.namedtuple(
            f'{source_db.capitalize()}To{target_db.capitalize()}Interaction',
            [
                f'{source_db}_name',
                f'{target_db.capitalize()}Entries',
                'chebi_id'
            ]
        )

    else:
        Interaction = collections.namedtuple(
            f'{source_db.capitalize()}To{target_db.capitalize()}Interaction',
            [
                f'{source_db}_name',
                f'{target_db.capitalize()}Entries',
            ]
        )

    source = db_list[0]
    target = db_list[1]

    source_url = source_db if source_db != 'gene' else org
    target_url = target_db if target_db != 'gene' else org

    entries = _kegg_link(source_url, target_url)
    interactions = dict()

    for entry in entries:

        source_id = source.handle(entry[0])

        try:
            source_name = source.get(source_id)
        except:
            source_name = None

        target_id = target.handle(entry[1])

        try:
            target_name = target.get(target_id)
        except:
            target_name = None

        if target_db == 'gene':

            ncbi_gene_id = kegg_to_ncbi.get(target_id)
            uniprot_ids = kegg_to_uniprot.get(target_id)

            if isinstance(ncbi_gene_id, list):
                ncbi_gene_id = tuple(ncbi_gene_id)

            if isinstance(uniprot_ids, list):
                uniprot_ids = tuple(uniprot_ids)

            target_db_entry = TargetDbEntry(
                target_id,
                target_name,
                ncbi_gene_id,
                uniprot_ids
            )
        
        elif target_db == 'drug':

            chebi_id = kegg_to_chebi.get(target_id)

            if isinstance(chebi_id, list):
                chebi_id = tuple(chebi_id)

            target_db_entry = TargetDbEntry(
                target_id,
                target_name,
                chebi_id
            )
        
        else:

            target_db_entry = TargetDbEntry(
                target_id,
                target_name,
            )

        try:
            interactions[source_id][f'{target_db}_entries'].append(target_db_entry)

        except KeyError:
            interactions[source_id] = dict()
            interactions[source_id][f'{source_db}_name'] = source_name
            interactions[source_id][f'{target_db}_entries'] = [target_db_entry]

            if source_db == 'gene':
                interactions[source_id]['ncbi_gene_id'] = kegg_to_ncbi.get(source_id)
                interactions[source_id]['uniprot_ids'] = kegg_to_uniprot.get(source_id)
            
            elif source_db == 'drug':
                interactions[source_id]['chebi_id'] = kegg_to_chebi.get(source_id)

    for key, value in interactions.items():

        if source_db == 'gene':
            interaction = Interaction (
                value[f'{source_db}_name'],
                tuple(value[f'{target_db}_entries']),
                value['ncbi_gene_id'],
                value['uniprot_ids']
            )
            
        elif source_db == 'drug':
            interaction = Interaction (
                value[f'{source_db}_name'],
                tuple(value[f'{target_db}_entries']),
                value['chebi_id'],
            )

        else:
            interaction = Interaction (
                value[f'{source_db}_name'],
                tuple(value[f'{target_db}_entries'])
            )

        interactions[key] = interaction

    if org != None:
        organism = _Organism()
        org_id, org_name = organism.get(org)
        interactions['org_id'] = org_id
        interactions['org_name'] = org_name

    return interactions


class _KeggDatabase(ABC):

    _data = None


    @abstractmethod
    def __init__(self):
        pass


    @abstractmethod
    def handle(self):
        pass


    @abstractmethod
    def download_data(self):
        pass


    def get(self, index):
        return self._data[index]


    def get_data(self):
        return self._data

    
class _Organism(_KeggDatabase):

    def __init__(self):
        self.download_data()


    def download_data(self):
        entries = _kegg_list("genome")

        data = {}

        for genome_id, description in entries:
            # description = "hsa; Homo sapiens (human)"
            org, org_name = description.split(";", 1)

            data[self.handle(org.strip())] = [
                genome_id,
                org_name.strip(),
            ]

        self._data = data

    def handle(self, org):
        return org


class _Gene(_KeggDatabase):

    def __init__(self, org):
        self.download_data(org)


    def download_data(self, org):

        entries = _kegg_list(org)
        
        gene_slice = [row[0] for row in entries]

        name_slice = [row[-1] for row in entries]
        name_slice = [name.split(';')[-1] for name in name_slice]
        name_slice = [name.strip(' ') for name in name_slice]

        entries = zip(gene_slice, name_slice)
        self._data = {self.handle(gene) : gene_name for (gene, gene_name) in entries}


    def handle(self, gene):

        return gene


class _Pathway(_KeggDatabase):

    def __init__(self, org=None):
        self.download_data()

    def download_data(self, org=None):

        if org != None:
        
            entries = _kegg_list('pathway', org)

        else:

            entries = _kegg_list('pathway')

        self._data = {self.handle(pathway) : pathway_name for (pathway, pathway_name) in entries}

    def handle(self, pathway):

        pathway_re = re.compile(r'\d+')
        pathway_id = pathway_re.search(pathway)

        return 'map' + pathway_id.group()


class _SplitDatabase(_KeggDatabase):
    def __init__(self, entry_url):
        self.download_data(entry_url)


    def download_data(self, entry_url):
        entries = _kegg_list(entry_url)
        self._data = {self.handle(entry) : entry_name for (entry, entry_name) in entries}


    def handle(self, entry):
        return entry.split(':')[1] if ":" in entry else entry


class _Disease(_SplitDatabase):

    def __init__(self):
        super().__init__('disease')


class _Drug(_SplitDatabase):

    def __init__(self):
        super().__init__('drug')

        
class _Compound(_SplitDatabase):

    def __init__(self):
        super().__init__('compound')


class _ConversionTable:

    def __init__(self):
        self._table = {}

    @abstractmethod
    def download_table(self):
        pass

    def get(self, index):
        return self._table.get(index)

    def get_table(self):
        return self._table


class _OrgTable(_ConversionTable):

    def __init__(self, org=None):
        super().__init__()

        if org is not None:
            self.download_table(org)


class _KeggToNcbi(_OrgTable):
    
    def download_table(self, org):
        table = _kegg_conv(org, 'ncbi-geneid', target_split=True)
        self._table.update(table)


class _NcbiToKegg(_OrgTable):

    def download_table(self, org):
        table = _kegg_conv('ncbi-geneid', org, source_split=True)
        self._table.update(table)


class _KeggToUniprot(_OrgTable):

    def download_table(self, org):
        table = _kegg_conv(org, 'uniprot', target_split=True)
        self._table.update(table)


class _UniprotToKegg(_OrgTable):

    def download_table(self, org):
        table = _kegg_conv('uniprot', org, source_split=True)
        self._table.update(table)


class _KeggToChebi(_ConversionTable):

    def download_table(self):
        table = _kegg_conv('drug', 'chebi', source_split=True, target_split=True)
        self._table = table


class _ChebiToKegg(_ConversionTable):

    def download_table(self):
        table = _kegg_conv('chebi', 'drug', source_split=True, target_split=True)
        self._table = table

@dataclass
class RequestResult:
    identifier: str
    endpoint: str
    success: bool
    retries: int
    status_code: Optional[int]
    error: Optional[str]
    data: list
